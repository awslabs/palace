# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import unittest

from trace_audit_contract import trace_groups


class TraceAuditContractTest(unittest.TestCase):
    def contract(self):
        return {"Version": 1, "SourceIndices": [3, 7, 11, 20, 25, 31],
                "ConstrainedIndices": [7, 25],
                "RequiredConductorInterfaces": [
                    {"Conductor": 1, "Interface": "MA"},
                    {"Conductor": 1, "Interface": "MS"}],
                "SemanticExcitationClasses": [
                    {"Name": "global", "Kind": "smooth-global", "Indices": [3]},
                    {"Name": "localized", "Kind": "localized-matching",
                     "Indices": [11, 20]},
                    {"Name": "conductor-1-ma", "Kind": "conductor-interface",
                     "Conductor": 1, "Interface": "MA", "Indices": [31]},
                    {"Name": "conductor-1-ms", "Kind": "conductor-interface",
                     "Conductor": 1, "Interface": "MS", "Indices": [20]},
                    {"Name": "fixed-combination", "Kind": "deterministic-combination",
                     "Indices": [3, 31], "Coefficients": [1.0, -0.5]}]}

    def test_arbitrary_source_count_partition_and_role_coverage(self):
        groups = trace_groups(self.contract())
        self.assertEqual(groups["free"], [3, 11, 20, 31])
        self.assertEqual(groups["constrained"], [7, 25])
        self.assertEqual(groups["semantic"]["localized"],
                         {"kind": "localized-matching", "indices": [11, 20]})
        self.assertFalse(groups["PhysicsQualified"])

    def test_duplicate_names_kinds_and_missing_semantic_coverage_fail(self):
        mutations = {}
        duplicate_name = self.contract()
        duplicate_name["SemanticExcitationClasses"][1]["Name"] = "global"
        mutations["duplicate-name"] = duplicate_name
        missing_kind = self.contract()
        del missing_kind["SemanticExcitationClasses"][0]["Kind"]
        mutations["missing-kind"] = missing_kind
        missing_role = self.contract()
        missing_role["SemanticExcitationClasses"] = [
            item for item in missing_role["SemanticExcitationClasses"]
            if item["Kind"] != "localized-matching"]
        mutations["missing-role"] = missing_role
        duplicate_pair = self.contract()
        duplicate_pair["SemanticExcitationClasses"].append({
            "Name": "duplicate-interface", "Kind": "conductor-interface",
            "Conductor": 1, "Interface": "MA", "Indices": [3]})
        mutations["duplicate-interface"] = duplicate_pair
        for name, contract in mutations.items():
            with self.subTest(name=name), self.assertRaises(ValueError):
                trace_groups(contract)

    def test_invalid_source_partition_interface_and_combination_fail(self):
        mutations = []
        empty = self.contract(); empty["SourceIndices"] = []
        mutations.append(empty)
        out_of_bank = self.contract(); out_of_bank["ConstrainedIndices"] = [99]
        mutations.append(out_of_bank)
        missing_interface = self.contract()
        missing_interface["SemanticExcitationClasses"] = [
            item for item in missing_interface["SemanticExcitationClasses"]
            if item.get("Interface") != "MS"]
        mutations.append(missing_interface)
        bad_coefficients = self.contract()
        next(item for item in bad_coefficients["SemanticExcitationClasses"]
             if item["Kind"] == "deterministic-combination")["Coefficients"] = [1]
        mutations.append(bad_coefficients)
        for contract in mutations:
            with self.assertRaises(ValueError):
                trace_groups(contract)


if __name__ == "__main__":
    unittest.main()
