# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import unittest

from trace_audit_contract import trace_groups


class TraceAuditContractTest(unittest.TestCase):
    def test_arbitrary_source_count_and_partition(self):
        contract = {"Version": 1, "SourceIndices": [3, 7, 11, 20, 25, 31],
                    "ConstrainedIndices": [7, 25],
                    "SemanticExcitationClasses": [
                        {"Name": "global", "Indices": [3]},
                        {"Name": "localized", "Indices": [11, 20]},
                        {"Name": "conductor-b", "Indices": [31]},
                    ]}
        groups = trace_groups(contract)
        self.assertEqual(groups["free"], [3, 11, 20, 31])
        self.assertEqual(groups["constrained"], [7, 25])
        self.assertEqual(groups["semantic"]["localized"], [11, 20])

    def test_empty_or_out_of_bank_contract_fails_closed(self):
        for contract in ({"Version": 1, "SourceIndices": [], "ConstrainedIndices": [],
                          "SemanticExcitationClasses": []},
                         {"Version": 1, "SourceIndices": [1, 2], "ConstrainedIndices": [3],
                          "SemanticExcitationClasses": [{"Name": "global", "Indices": [1]}]}):
            with self.assertRaises(ValueError):
                trace_groups(contract)


if __name__ == "__main__":
    unittest.main()
