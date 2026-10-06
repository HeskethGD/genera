# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Fixed

- Chart continuation now checks local sheet separation and full-step/half-step
  consistency, preventing the silent branch switch reported in issue #4.
  Ambiguous quadrature nodes are refined and checked against endpoint branch
  labels. Unresolved continuation raises `NoConvergence`.

## [0.1.0] - 2026-10-04

### Initial Release

This is an experimental research release providing arbitrary-precision computation with algebraic curves and Abelian functions.

**What's supported:**
- Algebraic curve representation with automatic period matrix computation
- Abel maps for integration on algebraic curves
- Kleinian sigma and P functions for Abelian function theory
- Abelian function evaluation at arbitrary precision using mpmath
- Arbitrary-precision arithmetic throughout

**Known limitations:**
- API may change in minor releases during the 0.x series
- Numerical stability edge cases under investigation
- Windows testing not yet comprehensive

**API stability promise:**
- 0.x releases are experimental
- Breaking changes may occur in 0.x minor versions with notice in changelog
- Will follow strict semantic versioning from 1.0.0 onwards
- Public API changes will be documented in this changelog

**Features:**
- Algebraic curve representation and period matrix computation
- Kleinian sigma functions with canonical hyperelliptic normalization
- Kleinian P functions (derivatives of log sigma)
- Abelian function evaluation at arbitrary argument vectors
- Integration with mpmath for arbitrary-precision numerics

[0.1.0]: https://github.com/HeskethGD/generapy/releases/tag/v0.1.0
