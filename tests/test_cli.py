import json

from ile_predict import cli


def test_cli_json_includes_status_and_policy(capsys):
    rc = cli.main(["verapamil", "--json", "--no-ci"])
    out = capsys.readouterr()
    assert rc == 0
    payload = json.loads(out.out)
    assert payload[0]["status"] == "bundled-lookup"
    assert payload[0]["confidence_policy"] in {"standard", "low"}


def test_cli_reports_low_confidence_for_domain_edge(capsys):
    cli._print_reliability_notes(
        [{"matched": "test-agent", "confidence_policy": "low"}]
    )
    out = capsys.readouterr()
    assert "[low confidence] test-agent" in out.err
