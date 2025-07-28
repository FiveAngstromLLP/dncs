#!/bin/bash

# Forcefield Accuracy Test Runner
# This script runs the comprehensive forcefield test suite for DNCS

set -e  # Exit on any error

echo "=================================================="
echo "DNCS Forcefield Accuracy Test Suite"
echo "=================================================="
echo

# Function to run tests with timing
run_test_suite() {
    local test_name=$1
    local test_file=$2
    
    echo "Running $test_name..."
    echo "--------------------------------------------------"
    
    start_time=$(date +%s)
    
    if cargo test --test "$test_file" -- --nocapture; then
        end_time=$(date +%s)
        duration=$((end_time - start_time))
        echo "✅ $test_name completed successfully in ${duration}s"
    else
        echo "❌ $test_name failed"
        exit 1
    fi
    
    echo
}

# Function to run specific test categories
run_test_category() {
    local category=$1
    local test_file=$2
    local pattern=$3
    
    echo "Running $category tests..."
    echo "--------------------------------------------------"
    
    if cargo test --test "$test_file" "$pattern" -- --nocapture; then
        echo "✅ $category tests completed successfully"
    else
        echo "❌ $category tests failed"
        exit 1
    fi
    
    echo
}

# Check if cargo is available
if ! command -v cargo &> /dev/null; then
    echo "❌ Error: cargo is not installed or not in PATH"
    exit 1
fi

# Check if we're in the right directory
if [ ! -f "Cargo.toml" ]; then
    echo "❌ Error: Not in a Rust project directory (Cargo.toml not found)"
    exit 1
fi

# Parse command line arguments
case "${1:-all}" in
    "all")
        echo "Running all forcefield accuracy tests..."
        echo
        
        # Run main accuracy tests
        run_test_suite "Main Forcefield Accuracy Tests" "forcefield_accuracy_tests"
        
        # Run component tests
        run_test_suite "Force Component Tests" "force_components_tests"
        
        echo "=================================================="
        echo "🎉 All tests completed successfully!"
        echo "=================================================="
        ;;
        
    "accuracy")
        run_test_suite "Forcefield Accuracy Tests" "forcefield_accuracy_tests"
        ;;
        
    "components")
        run_test_suite "Force Component Tests" "force_components_tests"
        ;;
        
    "energy")
        run_test_category "Energy Calculation" "forcefield_accuracy_tests" "test_energy"
        ;;
        
    "forces")
        run_test_category "Force Components" "force_components_tests" "test_harmonic"
        run_test_category "Non-bonded Forces" "force_components_tests" "test_lennard_jones"
        run_test_category "Electrostatic Forces" "force_components_tests" "test_electrostatic"
        ;;
        
    "performance")
        run_test_category "Performance" "forcefield_accuracy_tests" "test_performance"
        ;;
        
    "integration")
        run_test_category "Integration" "forcefield_accuracy_tests" "test_integration"
        ;;
        
    "quick")
        echo "Running quick test subset..."
        echo
        
        # Run a subset of critical tests
        run_test_category "Basic Energy Tests" "forcefield_accuracy_tests" "test_energy_calculation_consistency"
        run_test_category "Parameter Validation" "force_components_tests" "test_forcefield_parameters"
        run_test_category "Integration Test" "forcefield_accuracy_tests" "test_integration"
        
        echo "✅ Quick tests completed successfully"
        ;;
        
    "help"|"-h"|"--help")
        echo "Usage: $0 [test_category]"
        echo
        echo "Test categories:"
        echo "  all         - Run all tests (default)"
        echo "  accuracy    - Run main accuracy tests only"
        echo "  components  - Run force component tests only"
        echo "  energy      - Run energy calculation tests"
        echo "  forces      - Run force component tests"
        echo "  performance - Run performance tests"
        echo "  integration - Run integration tests"
        echo "  quick       - Run a quick subset of critical tests"
        echo "  help        - Show this help message"
        echo
        echo "Examples:"
        echo "  $0              # Run all tests"
        echo "  $0 accuracy     # Run only accuracy tests"
        echo "  $0 quick        # Run quick test subset"
        echo
        exit 0
        ;;
        
    *)
        echo "❌ Error: Unknown test category '$1'"
        echo "Use '$0 help' to see available options"
        exit 1
        ;;
esac