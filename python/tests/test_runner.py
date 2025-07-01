#!/usr/bin/env python3
"""
Test runner script for Palace Python package.

This script can be used to run tests independently of pytest command line.
"""

import argparse
import os
import subprocess
import sys


def run_unit_tests():
    """Run unit tests only."""
    print("Running unit tests...")
    cmd = [sys.executable, "-m", "pytest", "tests/unit/", "-v"]
    return subprocess.run(cmd).returncode


def run_integration_tests():
    """Run integration tests only."""
    print("Running integration tests...")
    cmd = [sys.executable, "-m", "pytest", "tests/integration/", "-v"]
    return subprocess.run(cmd).returncode


def run_all_tests():
    """Run all tests."""
    print("Running all tests...")
    cmd = [sys.executable, "-m", "pytest", "tests/", "-v"]
    return subprocess.run(cmd).returncode


def run_example_verification():
    """Run a quick verification that examples work."""
    print("Verifying examples can be imported and basic functions work...")

    examples_dir = os.path.join(os.path.dirname(__file__), "..", "examples")
    examples_dir = os.path.abspath(examples_dir)

    if not os.path.exists(examples_dir):
        print(f"Examples directory not found: {examples_dir}")
        return 1

    # Test basic_usage.py
    test_script = f"""
import sys
sys.path.insert(0, r'{examples_dir}')
sys.path.insert(0, r'{os.path.dirname(examples_dir)}')

try:
    import basic_usage
    print("✓ basic_usage.py imports successfully")

    # Test config creation (should work without Palace)
    basic_usage.example_3_create_config()
    print("✓ basic_usage config creation works")

except Exception as e:
    print(f"✗ basic_usage.py error: {{e}}")
    sys.exit(1)

try:
    import eigenmode_analysis
    print("✓ eigenmode_analysis.py imports successfully")

    # Test config creation
    config = eigenmode_analysis.create_cavity_eigenmode_config("test.msh")
    assert config['Problem']['Type'] == 'Eigenmode'
    print("✓ eigenmode_analysis config creation works")

except Exception as e:
    print(f"✗ eigenmode_analysis.py error: {{e}}")
    sys.exit(1)

try:
    import frequency_domain_simulation
    print("✓ frequency_domain_simulation.py imports successfully")

    # Test config creation
    config = frequency_domain_simulation.create_two_port_network_config("test.msh")
    assert config['Problem']['Type'] == 'Driven'
    print("✓ frequency_domain_simulation config creation works")

except Exception as e:
    print(f"✗ frequency_domain_simulation.py error: {{e}}")
    sys.exit(1)

try:
    import time_domain_simulation
    print("✓ time_domain_simulation.py imports successfully")

    # Test config creation
    config = time_domain_simulation.create_transient_config("test.msh")
    assert config['Problem']['Type'] == 'Transient'
    print("✓ time_domain_simulation config creation works")

except Exception as e:
    print(f"✗ time_domain_simulation.py error: {{e}}")
    sys.exit(1)

print("All example verifications passed!")
"""

    return subprocess.run([sys.executable, "-c", test_script]).returncode


def main():
    """Main test runner."""
    parser = argparse.ArgumentParser(description="Palace Python Test Runner")
    parser.add_argument("--unit", action="store_true", help="Run unit tests only")
    parser.add_argument(
        "--integration", action="store_true", help="Run integration tests only"
    )
    parser.add_argument(
        "--examples", action="store_true", help="Run example verification only"
    )
    parser.add_argument("--all", action="store_true", help="Run all tests (default)")

    args = parser.parse_args()

    # Change to the directory containing this script
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    os.chdir("..")  # Go to python/ directory

    exit_code = 0

    if args.unit:
        exit_code = run_unit_tests()
    elif args.integration:
        exit_code = run_integration_tests()
    elif args.examples:
        exit_code = run_example_verification()
    else:
        # Run all by default
        print("Running comprehensive test suite...")

        # First run example verification (fastest)
        print("\n" + "=" * 60)
        print("1. Example Verification")
        print("=" * 60)
        exit_code = run_example_verification()

        if exit_code == 0:
            # Then run unit tests
            print("\n" + "=" * 60)
            print("2. Unit Tests")
            print("=" * 60)
            exit_code = run_unit_tests()

        if exit_code == 0:
            # Finally run integration tests
            print("\n" + "=" * 60)
            print("3. Integration Tests")
            print("=" * 60)
            exit_code = run_integration_tests()

        if exit_code == 0:
            print("\n" + "=" * 60)
            print("🎉 ALL TESTS PASSED!")
            print("=" * 60)
        else:
            print("\n" + "=" * 60)
            print("❌ SOME TESTS FAILED!")
            print("=" * 60)

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
