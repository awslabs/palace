#!/usr/bin/env python3

import importlib.util
import tempfile
import unittest
from pathlib import Path


MODULE_DIRECTORY = Path(__file__).parent
SCALAR_SPEC = importlib.util.spec_from_file_location(
    "local_edge_correction", MODULE_DIRECTORY / "local_edge_correction.py"
)
SCALAR = importlib.util.module_from_spec(SCALAR_SPEC)
SCALAR_SPEC.loader.exec_module(SCALAR)

import sys

sys.modules["local_edge_correction"] = SCALAR
POLARIZED_SPEC = importlib.util.spec_from_file_location(
    "local_edge_polarized_correction",
    MODULE_DIRECTORY / "local_edge_polarized_correction.py",
)
POLARIZED = importlib.util.module_from_spec(POLARIZED_SPEC)
POLARIZED_SPEC.loader.exec_module(POLARIZED)


class LocalEdgePolarizedCorrectionTest(unittest.TestCase):
    @staticmethod
    def write_case(directory, *, fabricated):
        total_n = 9.0 if fabricated else 10.0
        total_t = 17.0 if fabricated else 20.0
        inside_n = 5.0 if fabricated else 6.0
        inside_t = 9.0 if fabricated else 12.0
        total = total_n + total_t
        (directory / "surface-Q.csv").write_text(
            "i,p_surf[1],Q_surf[1]\n"
            f"1.0,{total},100.0\n"
        )
        header = [
            "i",
            "exc",
            "interface",
            "edge",
            "R (m)",
            "p_total",
            "p_in",
            "p_ann",
            "p_bulk_ann",
            "p_total_n",
            "p_total_t",
            "p_in_n",
            "p_in_t",
            "p_ann_n",
            "p_ann_t",
            *(f"p_bulk_{channel}_ann" for channel in POLARIZED.CHANNELS),
        ]
        rows = []
        for radius in (2.0e-7, 3.0e-7, 5.0e-7):
            ann_n = 1.0 + 1.0e6 * radius
            ann_t = 2.0 + 2.0e6 * radius
            channels = {
                "top_n": 2.0 * radius,
                "top_m": 3.0 * radius,
                "top_l": 1.0 * radius,
                "bottom_n": 4.0 * radius,
                "bottom_m": 0.0,
                "bottom_l": 0.0,
            }
            row = [
                1.0,
                0,
                1,
                1,
                radius,
                total,
                inside_n + inside_t,
                ann_n + ann_t,
                sum(channels.values()),
                total_n,
                total_t,
                inside_n,
                inside_t,
                ann_n,
                ann_t,
                *(channels[channel] for channel in POLARIZED.CHANNELS),
            ]
            rows.append(",".join(str(value) for value in row))
        (directory / "surface-Q-edge-local.csv").write_text(
            ",".join(header) + "\n" + "\n".join(rows) + "\n"
        )
        (directory / "case_resolved.json").write_text(
            '{"Model": {"L0": 1.0}, "Solver": {"Order": 2}, '
            '"Boundaries": {"Postprocessing": {"Dielectric": ['
            '{"Index": 1, "Type": "SA", "Thickness": 2e-9, '
            '"Permittivity": 4.0, "EdgeFrameNormal": [0, 1, 0]}]}}}'
        )

    def test_descriptor_mapping_is_process_side_aware(self):
        self.assertEqual(POLARIZED.DESCRIPTOR_CHANNELS[("MA", "n")], ("top_n",))
        self.assertEqual(
            POLARIZED.DESCRIPTOR_CHANNELS[("MS", "n")], ("bottom_n",)
        )
        self.assertEqual(POLARIZED.DESCRIPTOR_CHANNELS[("SA", "t")], ("top_m",))

    def test_polarized_calibration_and_application_reconstruct_reference(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            thin = root / "thin"
            fabricated = root / "fabricated"
            thin.mkdir()
            fabricated.mkdir()
            self.write_case(thin, fabricated=False)
            self.write_case(fabricated, fabricated=True)

            calibration_path = root / "calibration.csv"
            calibration = POLARIZED.calibrate(
                thin, fabricated, calibration_path, fit_count=3
            )
            by_component = {
                row["surface_component"]: next(
                    candidate
                    for candidate in calibration
                    if candidate["surface_component"] == row["surface_component"]
                )
                for row in calibration
            }
            self.assertAlmostEqual(
                by_component["n"]["descriptor_amplitude_per_m"], 2.0
            )
            self.assertAlmostEqual(
                by_component["t"]["descriptor_amplitude_per_m"], 3.0
            )

            rows = POLARIZED.apply(
                calibration_path,
                thin,
                root / "corrected.csv",
                [],
                [0.2, 0.3, 0.5],
                3,
                0.1,
                fabricated,
            )
            summary = next(row for row in rows if row["method"] == "window-mean")
            self.assertAlmostEqual(summary["p_corrected_n"], 9.0)
            self.assertAlmostEqual(summary["p_corrected_t"], 17.0)
            self.assertAlmostEqual(summary["p_corrected"], 26.0)
            self.assertAlmostEqual(summary["relative_error"], 0.0)
            self.assertAlmostEqual(
                summary["unmodeled_longitudinal_fraction"], 0.25
            )
            self.assertAlmostEqual(summary["model_confidence"], 0.75)

    def test_reader_rejects_nonconservative_bulk_channels(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            self.write_case(directory, fabricated=False)
            path = directory / "surface-Q-edge-local.csv"
            text = path.read_text().replace("2e-06,10.0,20.0", "3e-06,10.0,20.0", 1)
            path.write_text(text)
            with self.assertRaisesRegex(RuntimeError, "do not conserve p_bulk_ann"):
                POLARIZED.read_local_edges(directory)


if __name__ == "__main__":
    unittest.main()
