#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Every example configuration with an AMR size cap states the predictive size check.

Model.Refinement.MaxSizePredicted is false by default (upstream Palace's rule): AMR then
solves one more level above MaxSize (measured on the cpw2d calibration examples: +64-78 %
unknowns and 1.8-2.5x the peak memory of the last level; review r2152 C s12 item 1). So an
example under examples/ whose Model.Refinement has MaxSize > 0 must say which rule it uses:
MaxSizePredicted is present (true or false) beside MaxSize. Palace's configuration files
may hold // and /* */ comments and trailing commas (the reader ignores them); they are
removed here before the JSON is parsed.
"""

import json
import unittest
from pathlib import Path

EXAMPLES = Path(__file__).resolve().parent


def _strip(text, comments, trailing_commas):
    """Remove comments or trailing commas outside JSON strings."""
    out = []
    i, n, in_string = 0, len(text), False
    while i < n:
        c = text[i]
        if in_string:
            out.append(c)
            if c == "\\" and i + 1 < n:
                out.append(text[i + 1])
                i += 2
                continue
            in_string = c != '"'
            i += 1
            continue
        if c == '"':
            in_string = True
        elif comments and text.startswith("//", i):
            j = text.find("\n", i)
            i = n if j < 0 else j
            continue
        elif comments and text.startswith("/*", i):
            j = text.find("*/", i + 2)
            i = n if j < 0 else j + 2
            continue
        elif trailing_commas and c == ",":
            j = i + 1
            while j < n and text[j] in " \t\r\n":
                j += 1
            if j < n and text[j] in "}]":
                i += 1
                continue
        out.append(c)
        i += 1
    return "".join(out)


def load_config(path):
    text = _strip(path.read_text(), comments=True, trailing_commas=False)
    return json.loads(_strip(text, comments=False, trailing_commas=True))


def capped_refinements(root):
    """(path, Model.Refinement) of every example config with Refinement.MaxSize > 0."""
    for path in sorted(root.rglob("*.json")):
        config = load_config(path)
        model = config.get("Model") if isinstance(config, dict) else None
        refinement = model.get("Refinement") if isinstance(model, dict) else None
        if isinstance(refinement, dict) and refinement.get("MaxSize", 0) > 0:
            yield path, refinement


class ExampleMaxSizePredictedTest(unittest.TestCase):
    def test_every_capped_example_states_max_size_predicted(self):
        capped = list(capped_refinements(EXAMPLES))
        # The six configurations of the review's list are capped; a parse that found none
        # would make the check vacuous.
        self.assertGreaterEqual(len(capped), 6)
        missing = [str(path.relative_to(EXAMPLES)) for path, refinement in capped
                   if "MaxSizePredicted" not in refinement]
        self.assertEqual(missing, [], "Model.Refinement has MaxSize > 0 without MaxSizePredicted")

    def test_the_reader_handles_comments_and_trailing_commas(self):
        text = '{"a": "x // y", /* c */ "b": [1, 2,], // d\n "c": {"e": 1,},}'
        self.assertEqual(json.loads(_strip(_strip(text, True, False), False, True)),
                         {"a": "x // y", "b": [1, 2], "c": {"e": 1}})

    def test_a_capped_config_without_the_key_is_found(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, "capped.json").write_text(
                '{"Model": {"Refinement": {"MaxIts": 2, "MaxSize": 10, // cap\n}}}')
            Path(tmp, "uncapped.json").write_text('{"Model": {"Refinement": {"MaxIts": 2}}}')
            capped = list(capped_refinements(Path(tmp)))
        self.assertEqual([path.name for path, _ in capped], ["capped.json"])
        self.assertNotIn("MaxSizePredicted", capped[0][1])


if __name__ == "__main__":
    unittest.main()
