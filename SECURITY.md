# Security Policy

## Scope

This repository contains local numerical-analysis examples. It is not a network service and does not intentionally process authentication data, personal data, or remote commands.

## Safe-use boundary

Run the tools only on input files you trust. In particular, STEP parsing depends on third-party CAD and geometry libraries. A malicious CAD file could exploit a vulnerability in a dependency even though the repository itself does not execute commands or make network requests.

For untrusted geometry:

1. use a disposable virtual environment;
2. run without administrator privileges;
3. work in a disposable directory;
4. keep the Python dependencies patched;
5. do not place secrets or sensitive files in the working directory;
6. review generated files before sharing them.

The scripts should not be granted access to production systems, credentials, private project directories, or connected machinery.

## Reporting a vulnerability

Do not publish an exploitable vulnerability in a public issue. Contact the maintainer privately through the email address listed on the maintainer's GitHub profile, including:

- affected file and version/commit;
- reproduction steps;
- impact;
- a minimal proof of concept that does not contain secrets.

Security fixes will be disclosed after a fix or mitigation is available.
