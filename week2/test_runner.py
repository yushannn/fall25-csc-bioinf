#!/usr/bin/env python3
"""
Test runner script for week2 CI validation.

This script runs all tests and provides comprehensive error reporting
suitable for CI environments.
"""

import sys
import os
import traceback
import subprocess
from typing import List, Tuple

def run_test_script(script_name: str, description: str) -> Tuple[bool, str]:
    """Run a test script and capture its output.
    
    Args:
        script_name: Name of the Python script to run
        description: Human-readable description of the test
        
    Returns:
        Tuple of (success, output_message)
    """
    try:
        print(f"\n{'='*60}")
        print(f"Running: {description}")
        print(f"Script: {script_name}")
        print('='*60)
        
        # Run the script
        result = subprocess.run(
            [sys.executable, script_name],
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )
        
        # Print the output
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr, file=sys.stderr)
            
        if result.returncode == 0:
            print(f"✅ {description} PASSED")
            return True, f"{description} passed successfully"
        else:
            print(f"❌ {description} FAILED (exit code: {result.returncode})")
            return False, f"{description} failed with exit code {result.returncode}"
            
    except subprocess.TimeoutExpired:
        error_msg = f"{description} timed out after 5 minutes"
        print(f"⏰ {error_msg}")
        return False, error_msg
    except Exception as e:
        error_msg = f"{description} failed with exception: {str(e)}"
        print(f"💥 {error_msg}")
        traceback.print_exc()
        return False, error_msg


def main():
    """Run all tests and report results."""
    print("Week 2 Bio.motifs Test Runner")
    print("="*60)
    
    # Define test suite
    tests = [
        ("test.py", "BioPython Compatibility Tests"),
        ("test_comprehensive.py", "Comprehensive Test Suite"),
        ("benchmark_performance.py", "Performance Benchmarks"),
    ]
    
    results = []
    total_tests = len(tests)
    passed_tests = 0
    
    # Run each test
    for script, description in tests:
        if not os.path.exists(script):
            print(f"⚠️  Warning: {script} not found, skipping {description}")
            results.append((False, f"{description} - script not found"))
            continue
            
        success, message = run_test_script(script, description)
        results.append((success, message))
        if success:
            passed_tests += 1
    
    # Final report
    print("\n" + "="*60)
    print("FINAL TEST RESULTS")
    print("="*60)
    
    for i, (success, message) in enumerate(results, 1):
        status = "✅ PASS" if success else "❌ FAIL"
        print(f"{i}. {status}: {message}")
    
    print(f"\nSummary: {passed_tests}/{total_tests} tests passed")
    
    if passed_tests == total_tests:
        print("\n🎉 All tests passed! Week 2 implementation is working correctly.")
        sys.exit(0)
    else:
        print(f"\n💥 {total_tests - passed_tests} test(s) failed. Please check the output above.")
        sys.exit(1)


if __name__ == "__main__":
    main()