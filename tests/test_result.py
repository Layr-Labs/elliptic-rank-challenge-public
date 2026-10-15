"""The `internal error: <Type>: <message>` detail is the one shape every entry point prints for a
non-VerifyError exception (GP_ERROR), so agents can parse it; one helper owns the format."""

from ec_verify.result import internal_error_detail


def test_internal_error_detail_names_the_type_and_truncates_the_message():
    prefix = "internal error: KeyError: "
    d = internal_error_detail(KeyError("x" * 500))
    assert d.startswith(prefix) and len(d) == len(prefix) + 200
    assert internal_error_detail(ValueError("short")) == "internal error: ValueError: short"
