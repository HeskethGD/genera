"""Sphinx directive for running a Genera example during documentation builds."""

from pathlib import Path
import os
import subprocess
import sys

from docutils import nodes
from docutils.parsers.rst import Directive, directives


class GeneraExampleDirective(Directive):
    """Run an example module and render its captured output."""

    required_arguments = 1
    optional_arguments = 0
    option_spec = {
        "caption": directives.unchanged,
        "timeout": directives.nonnegative_int,
    }

    def run(self):
        module = self.arguments[0]
        root = Path(__file__).resolve().parents[1]
        environment = os.environ.copy()
        paths = [str(root / "src"), str(root)]
        if environment.get("PYTHONPATH"):
            paths.append(environment["PYTHONPATH"])
        environment["PYTHONPATH"] = os.pathsep.join(paths)
        timeout = self.options.get("timeout", 120)
        try:
            completed = subprocess.run(
                [sys.executable, "-m", module],
                cwd=root,
                env=environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=timeout,
                check=False,
            )
        except subprocess.TimeoutExpired as exc:
            raise self.error(
                f"Genera example {module!r} exceeded the {timeout}-second "
                "documentation timeout"
            ) from exc

        output = completed.stdout or "(example produced no output)\n"
        if completed.returncode:
            raise self.error(
                f"Genera example {module!r} failed with exit code "
                f"{completed.returncode}:\n{output}"
            )

        result = []
        caption = self.options.get("caption")
        if caption:
            result.append(nodes.paragraph(text=caption))
        block = nodes.literal_block(output, output)
        block["language"] = "text"
        result.append(block)
        return result


def setup(app):
    app.add_directive("genera-example", GeneraExampleDirective)
    return {
        "version": "1",
        "parallel_read_safe": False,
    }
