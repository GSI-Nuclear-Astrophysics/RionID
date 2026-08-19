# Security Policy

## Reporting a vulnerability

If you discover a security vulnerability in RionID, please report it
privately by emailing D.FreireFernandez@gsi.de rather than opening a
public issue. Include a description of the vulnerability, steps to
reproduce it, and its potential impact.

You should expect an initial response within 5 business days. This is a
small, academically-maintained research-software project without a formal
security team — response times may vary.

## Scope

RionID is a desktop GUI/CLI application for offline scientific data
analysis. It does not run as a network service, does not handle
authentication or user accounts, and does not process untrusted network
input by design. Its main external-data-handling surface is:

- Reading user-supplied spectrum/candidate files (`.npz`/`.csv`/`.lpp`/
  binary formats) — a malformed file could in principle cause a crash or
  resource-exhaustion condition, not remote code execution.
- Downloading the public AME2020 mass table from IAEA
  (`www-nds.iaea.org`) over HTTPS on first use.

Dependency vulnerabilities are monitored via `pip-audit` in CI and
Dependabot's automated update PRs (`.github/dependabot.yml`).

## Supported versions

Only the latest released version is supported with security fixes.
