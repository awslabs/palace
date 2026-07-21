#!/usr/bin/env python3

import importlib.util
import math
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("local_edge_correction.py")
SPEC = importlib.util.spec_from_file_location("local_edge_correction", MODULE_PATH)
CORRECTION = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CORRECTION)


class LocalEdgeCorrectionTest(unittest.TestCase):
    @staticmethod
    def write_case(directory, inside, outside, *, order=2):
        total = inside + outside
        (directory / "surface-Q.csv").write_text(
            "i,p_surf[1],Q_surf[1]\n"
            f"1.0,{total},100.0\n"
        )
        rows = []
        for radius in (2.0e-7, 3.0e-7, 5.0e-7):
            surface_annulus = 1.0 + 1.0e6 * radius
            bulk_annulus = 2.0 * radius
            rows.append(
                f"1.0,0,1,1,{radius},{total},{inside},"
                f"{surface_annulus},{bulk_annulus}"
            )
        (directory / "surface-Q-edge-local.csv").write_text(
            "i,exc,interface,edge,R (m),p_total,p_in,p_ann,p_bulk_ann\n"
            + "\n".join(rows)
            + "\n"
        )
        (directory / "case_resolved.json").write_text(
            '{"Model": {"L0": 1.0}, "Solver": {"Order": '
            f"{order}"
            '}, "Boundaries": {"Postprocessing": {"Dielectric": '
            '[{"Index": 1, "Type": "SA", "Thickness": 2e-9, '
            '"Permittivity": 4.0}]}}}'
        )

    @staticmethod
    def calibration_rows():
        return [
            {
                "interface": 1,
                "radius_um": radius_um,
                "inside_singular_coefficient_m": 1.0,
                "outside_delta_coefficient_m": 0.5,
                "inside_coefficient_amr_spread": 0.0,
                "outside_delta_amr_spread": 0.0,
            }
            for radius_um in (0.2, 0.3, 0.5)
        ]

    @staticmethod
    def window_candidate(radius, corrected):
        return {
            "method": "radius",
            "radius_um": radius,
            "p_corrected": corrected,
            "Q_corrected": 100.0,
            "p_reference": 2.0,
            "p_out": 0.5,
            "p_smooth_in": 0.25,
            "p_modeled_singular": 1.0,
            "p_outside_delta": 0.25,
            "singular_fraction": 1.0,
            "singular_identifiability": 1.0,
            "fit_confidence": 1.0,
            "modeled_fraction": 0.625,
            "model_confidence": 1.0,
            "fit_relative_rms": 0.0,
            "bulk_fit_relative_rms": 0.0,
            "surface_fit_relative_rms": 0.0,
            "sampled_edges": 1,
            "sampled_edge_fraction": 1.0,
            "calibration_inside_amr_spread": 0.0,
            "calibration_outside_amr_spread": 0.0,
        }

    def test_canonical_state_normalizes_numeric_spellings(self):
        self.assertEqual(CORRECTION.canonical_state("1.0"), 1.0)
        self.assertEqual(CORRECTION.canonical_state("1.00e+00"), 1.0)
        self.assertEqual(CORRECTION.canonical_state("mode"), "mode")

    def test_automatic_segments_pool_by_physical_chain(self):
        radius = 1.0
        values = {}
        for edge, length, vertex_types, participation in (
            (1, 2.0, (1, 0), 3.0),
            (2, 3.0, (0, 2), 4.0),
        ):
            values[(1.0, 0, 1, edge, radius)] = {
                "p_total": participation,
                "p_in": 0.5 * participation,
                "p_ann": 0.25 * participation,
                "p_bulk_ann": radius * participation,
                "automatic": 1,
                "component": 4,
                "chain": 7,
                "vertex_types": vertex_types,
                "edge_length": length,
            }

        groups = CORRECTION.group_local_edges(values)
        self.assertEqual(len(groups), 1)
        edge_values = {key[3]: value for key, value in groups.items()}
        pooled = next(iter(edge_values.values()))[radius]
        self.assertEqual(pooled["p_total"], 7.0)
        self.assertEqual(pooled["segment_count"], 2)
        self.assertEqual(pooled["edge_length"], 5.0)
        self.assertEqual(pooled["nonregular_endpoint_count"], 2)

        topology = CORRECTION.summarize_edge_topology(edge_values, radius)
        self.assertEqual(topology["automatic_chain_count"], 1)
        self.assertEqual(topology["nonregular_chain_count"], 1)
        self.assertEqual(topology["segment_count"], 2)
        self.assertAlmostEqual(
            topology["unmodeled_vertex_length_fraction"], 0.4
        )
        self.assertAlmostEqual(
            topology["unmodeled_vertex_surface_fraction"], 17.0 / 42.0
        )
        self.assertAlmostEqual(
            topology["unmodeled_vertex_bulk_fraction"], 17.0 / 42.0
        )
        self.assertAlmostEqual(topology["unmodeled_vertex_fraction"], 17.0 / 42.0)

    def test_exact_vertex_energies_override_endpoint_length_estimate(self):
        radius = 1.0
        values = {
            (1.0, 0, 1, 1, radius): {
                "p_total": 120.0,
                "p_in": 100.0,
                "p_ann": 30.0,
                "p_bulk_ann": 50.0,
                "p_vertex_in": 7.0,
                "p_bulk_vertex_ann": 3.0,
                "automatic": 1,
                "component": 4,
                "chain": 7,
                "vertex_types": (1, 2),
                "edge_length": 10.0,
            }
        }

        groups = CORRECTION.group_local_edges(values)
        edge_values = {key[3]: value for key, value in groups.items()}
        pooled = next(iter(edge_values.values()))[radius]
        self.assertEqual(pooled["vertex_risk_p_inside"], 7.0)
        self.assertEqual(pooled["vertex_risk_p_bulk_ann"], 3.0)

        topology = CORRECTION.summarize_edge_topology(edge_values, radius)
        self.assertAlmostEqual(
            topology["unmodeled_vertex_length_fraction"], 0.2
        )
        self.assertAlmostEqual(
            topology["unmodeled_vertex_surface_fraction"], 0.07
        )
        self.assertAlmostEqual(
            topology["unmodeled_vertex_bulk_fraction"], 0.06
        )

    def test_singular_fit_recovers_intercept(self):
        amplitude = 2.5
        slope = 4.0e5
        values = {
            radius: {"p_bulk_ann": radius * (amplitude + slope * radius)}
            for radius in (2.0e-7, 3.0e-7, 5.0e-7)
        }
        fitted, residual, points = CORRECTION.fit_singular_amplitude(values, 3)
        self.assertAlmostEqual(fitted, amplitude)
        self.assertLess(residual, 1.0e-14)
        self.assertEqual(len(points), 3)

    def test_surface_annulus_fit_recovers_smooth_density(self):
        intercept = 1.5
        smooth_density = 4.0e5
        values = {
            radius: {"p_ann": intercept + smooth_density * radius}
            for radius in (2.0e-7, 3.0e-7, 5.0e-7)
        }
        fitted_intercept, fitted_density, residual, points = (
            CORRECTION.fit_surface_annulus(values, 3)
        )
        self.assertAlmostEqual(fitted_intercept, intercept)
        self.assertAlmostEqual(fitted_density, smooth_density)
        self.assertLess(residual, 1.0e-14)
        self.assertEqual(len(points), 3)

    def test_smooth_field_can_have_high_model_confidence(self):
        singular, confidence = CORRECTION.edge_model_confidence(
            0.0, 5.0, 0.0, 0.0, 0.1
        )
        self.assertEqual(singular, 0.0)
        self.assertEqual(confidence, 1.0)

    def test_failed_fit_has_zero_model_confidence(self):
        singular, confidence = CORRECTION.edge_model_confidence(
            2.0, 4.0, math.inf, 0.0, 0.1
        )
        self.assertEqual(singular, 0.5)
        self.assertEqual(confidence, 0.0)

    def test_model_confidence_is_bounded_when_fit_exceeds_observation(self):
        singular, confidence = CORRECTION.edge_model_confidence(
            10.0, 2.0, 0.1, 0.05, 0.1
        )
        self.assertEqual(singular, 1.0)
        self.assertEqual(confidence, 0.5)

    def test_additive_correction_reports_each_decomposed_term(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            self.write_case(directory, inside=6.0, outside=4.0)
            rows = CORRECTION.calculate_radius_rows(
                self.calibration_rows(),
                directory,
                [],
                [],
                3,
                0.1,
                None,
            )
            self.assertEqual(len(rows), 3)
            row = rows[0]
            self.assertAlmostEqual(row["p_out"], 4.0)
            self.assertAlmostEqual(row["p_smooth_in"], 0.2)
            self.assertAlmostEqual(row["p_modeled_singular"], 2.0)
            self.assertAlmostEqual(row["p_outside_delta"], 1.0)
            self.assertAlmostEqual(row["p_corrected"], 7.2)
            self.assertAlmostEqual(row["model_confidence"], 1.0)

    def test_correction_does_not_retain_raw_inside_energy(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = root / "first"
            second = root / "second"
            first.mkdir()
            second.mkdir()
            self.write_case(first, inside=6.0, outside=4.0)
            self.write_case(second, inside=60.0, outside=4.0)
            first_rows = CORRECTION.calculate_radius_rows(
                self.calibration_rows(), first, [], [], 3, 0.1, None
            )
            second_rows = CORRECTION.calculate_radius_rows(
                self.calibration_rows(), second, [], [], 3, 0.1, None
            )
            self.assertNotEqual(first_rows[0]["p_raw"], second_rows[0]["p_raw"])
            for first_row, second_row in zip(first_rows, second_rows):
                self.assertAlmostEqual(
                    first_row["p_corrected"], second_row["p_corrected"]
                )

    def test_correction_confidence_accounts_for_singular_identifiability(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            self.write_case(directory, inside=6.0, outside=4.0)
            rows = []
            for radius in (2.0e-7, 3.0e-7, 5.0e-7):
                surface_annulus = 1.0 + 1.0e6 * radius
                bulk_annulus = radius * (2.0 + 4.0e7 * radius)
                rows.append(
                    f"1.0,0,1,1,{radius},10.0,6.0,"
                    f"{surface_annulus},{bulk_annulus}"
                )
            (directory / "surface-Q-edge-local.csv").write_text(
                "i,exc,interface,edge,R (m),p_total,p_in,p_ann,p_bulk_ann\n"
                + "\n".join(rows)
                + "\n"
            )
            calibration = self.calibration_rows()
            for row in calibration:
                row["inside_singular_coefficient_m"] = 10.0
                row["outside_delta_coefficient_m"] = 0.0
            result = CORRECTION.calculate_radius_rows(
                calibration, directory, [], [], 3, 0.1, None
            )[0]
            self.assertAlmostEqual(result["fit_confidence"], 1.0)
            self.assertAlmostEqual(result["singular_fraction"], 0.2)
            self.assertLess(result["model_confidence"], 0.4)

    def test_calibration_reconstructs_fabricated_reference(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            thin = root / "thin"
            fabricated = root / "fabricated"
            thin.mkdir()
            fabricated.mkdir()
            self.write_case(thin, inside=6.0, outside=4.0)
            self.write_case(fabricated, inside=3.5, outside=4.5)
            calibration = CORRECTION.calibrate(
                thin, fabricated, None, 3, include_history=False
            )
            rows = CORRECTION.calculate_radius_rows(
                calibration, thin, [], [], 3, 0.1, fabricated
            )
            for row in rows:
                self.assertAlmostEqual(row["p_corrected"], 8.0)
                self.assertAlmostEqual(row["relative_error"], 0.0)

    def test_local_assigned_participation_must_conserve_interface_total(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "surface-Q.csv").write_text(
                "i,p_surf[1],Q_surf[1]\n"
                "1.0,2.0,100.0\n"
            )
            rows = [
                f"1.0,0,1,1,{radius},1.5,0.5,{1.0 + 1.0e6 * radius},{radius}"
                for radius in (2.0e-7, 3.0e-7, 5.0e-7)
            ]
            (directory / "surface-Q-edge-local.csv").write_text(
                "i,exc,interface,edge,R (m),p_total,p_in,p_ann,p_bulk_ann\n"
                + "\n".join(rows)
                + "\n"
            )
            with self.assertRaisesRegex(
                RuntimeError, "Local assigned participations sum"
            ):
                CORRECTION.calculate_radius_rows(
                    self.calibration_rows(),
                    directory,
                    [],
                    [],
                    3,
                    0.1,
                    None,
                )

    def test_window_selection_prefers_smallest_stable_window(self):
        small = [{"radius_um": 0.2}]
        medium = [{"radius_um": 0.3}]
        large = [{"radius_um": 0.5}]
        selected = CORRECTION.choose_window(
            [(small, 0.03), (medium, 0.015), (large, 0.005)], 0.02
        )
        self.assertIs(selected[0], medium)

    def test_window_selection_falls_back_to_least_varying(self):
        small = [{"radius_um": 0.2}]
        large = [{"radius_um": 0.5}]
        selected = CORRECTION.choose_window(
            [(small, 0.06), (large, 0.03)], 0.02
        )
        self.assertIs(selected[0], large)

    def test_window_summary_uses_mean_instead_of_zero_radius_extrapolation(self):
        candidates = [
            self.window_candidate(radius, corrected)
            for radius, corrected in ((0.2, 1.0), (0.3, 2.0), (0.5, 3.0))
        ]
        summary = CORRECTION.summarize_window(candidates)
        self.assertEqual(summary["method"], "window-mean")
        self.assertAlmostEqual(summary["p_corrected"], 2.0)
        self.assertAlmostEqual(summary["relative_uncertainty"], 0.5)
        self.assertAlmostEqual(summary["relative_error"], 0.0)

    def test_window_uncertainty_includes_low_fit_confidence(self):
        candidate = self.window_candidate(0.2, 1.0)
        candidate["p_reference"] = math.nan
        candidate["singular_fraction"] = 0.1
        candidate["model_confidence"] = 0.1
        candidate["sampled_edge_fraction"] = 0.5
        candidates = [candidate]
        summary = CORRECTION.summarize_window(candidates)
        self.assertAlmostEqual(summary["relative_uncertainty"], 0.9)

    def test_read_fem_order_rejects_conflicting_configs(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "first_resolved.json").write_text(
                '{"Solver": {"Order": 1}}'
            )
            (directory / "second_resolved.json").write_text(
                '{"Solver": {"Order": 2}}'
            )
            with self.assertRaisesRegex(RuntimeError, "Conflicting FEM orders"):
                CORRECTION.read_fem_order(directory)

    def test_interface_metadata_uses_physical_thickness(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "case_resolved.json").write_text(
                '{"Model": {"L0": 1e-6}, "Boundaries": {"Postprocessing": '
                '{"Dielectric": [{"Index": 2, "Type": "MS", "Thickness": 0.002, '
                '"Permittivity": 11.47}]}}}'
            )
            metadata = CORRECTION.read_interface_metadata(directory)
            self.assertEqual(metadata[2]["interface_type"], "MS")
            self.assertAlmostEqual(metadata[2]["interface_thickness_m"], 2.0e-9)
            self.assertEqual(metadata[2]["interface_permittivity"], 11.47)
            self.assertFalse(metadata[2]["flux_recovery"])

    def test_iteration_directory_uses_parent_resolved_config(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            iteration = directory / "iteration3"
            iteration.mkdir()
            (directory / "case_resolved.json").write_text(
                '{"Solver": {"Order": 2}}'
            )
            self.assertEqual(CORRECTION.read_fem_order(iteration), 2)

    def test_interface_metadata_mismatch_is_rejected(self):
        calibration = {
            1: {
                "interface_type": "SA",
                "interface_thickness_m": 2.0e-9,
                "interface_permittivity": 4.0,
                "flux_recovery": True,
            }
        }
        target = {
            4: {
                "interface_type": "MS",
                "interface_thickness_m": 2.0e-9,
                "interface_permittivity": 11.47,
                "flux_recovery": True,
            }
        }
        with self.assertRaisesRegex(RuntimeError, "does not match calibration"):
            CORRECTION.validate_interface_metadata(calibration, target, 4, 1)

    def test_flux_recovery_mismatch_is_rejected(self):
        calibration = {
            1: {
                "interface_type": "SA",
                "interface_thickness_m": 2.0e-9,
                "interface_permittivity": 4.0,
                "flux_recovery": True,
            }
        }
        target = {
            4: {
                "interface_type": "SA",
                "interface_thickness_m": 2.0e-9,
                "interface_permittivity": 4.0,
                "flux_recovery": False,
            }
        }
        with self.assertRaisesRegex(
            RuntimeError, "flux recovery False != True"
        ):
            CORRECTION.validate_interface_metadata(calibration, target, 4, 1)

    def test_application_rejects_order_mismatch_before_reading_results(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            calibration = directory / "calibration.csv"
            calibration.write_text(
                "interface,radius_um,fem_order,interface_type,"
                "interface_thickness_m,interface_permittivity,"
                "flux_recovery,"
                "inside_singular_coefficient_m,outside_delta_coefficient_m,"
                "fit_count,inside_coefficient_amr_spread,"
                "outside_delta_amr_spread\n"
                "1,0.2,2,SA,2e-9,4.0,true,1e-9,1e-9,3,0.0,0.0\n"
            )
            target = directory / "target"
            target.mkdir()
            (target / "target_resolved.json").write_text(
                '{"Solver": {"Order": 1}}'
            )
            with self.assertRaisesRegex(
                RuntimeError, "Calibration uses FEM order p=2"
            ):
                CORRECTION.apply(
                    calibration,
                    target,
                    directory / "output.csv",
                    [],
                    [],
                    None,
                    0.1,
                    None,
                )

    def test_application_rejects_fit_count_mismatch(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            calibration = directory / "calibration.csv"
            calibration.write_text(
                "interface,radius_um,fem_order,interface_type,"
                "interface_thickness_m,interface_permittivity,"
                "flux_recovery,"
                "inside_singular_coefficient_m,outside_delta_coefficient_m,"
                "fit_count,inside_coefficient_amr_spread,"
                "outside_delta_amr_spread\n"
                "1,0.2,2,SA,2e-9,4.0,true,1e-9,1e-9,3,0.0,0.0\n"
            )
            with self.assertRaisesRegex(RuntimeError, "uses 3 fit radii"):
                CORRECTION.apply(
                    calibration,
                    directory,
                    directory / "output.csv",
                    [],
                    [],
                    4,
                    0.1,
                    None,
                )


if __name__ == "__main__":
    unittest.main()
