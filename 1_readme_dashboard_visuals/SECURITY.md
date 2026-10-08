# Security and privacy policy

This repository processes public health data. The things we most want to hear about are **privacy problems**, not
only software vulnerabilities.

## What to report privately

* **Personal data in the repository or its history**: any row, file, screenshot or log that identifies a person or
  could help re-identify one, including per-person rows from a source whose terms forbid redistribution.
* **A public source that turns out to contain direct identifiers** (names, dates of birth, record numbers, contact
  details). Two such records were found during development and deliberately not downloaded (figshare 30127723 and
  Mendeley vzydrgtymt; `docs/LAB_DATASET_DISCOVERY.md`); if you find another, report it here and, where appropriate,
  to the repository that hosts it. Do not download it again to show us.
* **Credentials or signed URLs** committed by mistake (API keys, tokens, pre-signed cloud storage links).
* **Vulnerabilities** in the API, dashboard, MCP server or container setup (for example path traversal through a tool
  argument, the services binding to a public interface by default, or the container running with more privilege than
  `docker-compose.yml` intends).
* **A guardrail bypass**: a tool, API route or MCP call that joins a participant to a place or facility, returns a
  participant's rows through `trace_evidence`, or presents a generated number as data.

## How to report

* **Use GitHub's private vulnerability reporting** on this repository (Security tab, "Report a vulnerability") once it
  is enabled, or
* e-mail the maintainer at **SECURITY_CONTACT_EMAIL** (to be filled in before the repository is made public).

Please include what you found, where (file and commit, URL or API route), and how to reproduce it. **Do not include the
personal data itself**: point to its location instead. We aim to acknowledge a report within 7 days.

## What not to do

* **Never put personal health information in a public issue, pull request, discussion or commit** — not your own, not
  anyone else's. Replace real values with made-up ones. Issues that contain it will be deleted.
* Do not open a public issue for a privacy or security problem before it has been fixed.
* Do not attempt to re-identify participants in any dataset to prove a point. The NCHS Data User Agreement, the SDR
  terms of the Stanford releases and this project's guardrails all forbid it.

## Scope and supported versions

This is an alpha research prototype; only the latest commit on the default branch is supported. The services are
designed to run locally and read-only (`127.0.0.1` by default, `MEASURE_IT_OFFLINE=1` in the container); exposing them
on a network is at your own risk. Nothing in this repository is a medical device or gives medical advice.
