#!/usr/bin/env python3
"""
Basic math functions for testing purposes.
"""

def add(a, b):
    """Add two numbers."""
    return a + b

def subtract(a, b):
    """Subtract two numbers."""
    return a - b

def multiply(a, b):
    """Multiply two numbers."""
    return a * b

def divide(a, b):
    """Divide two numbers."""
    if b == 0:
        return "Error: Division by zero"
    return a / b

def power(a, b):
    """Raise a to the power of b."""
    return a ** b

def modulo(a, b):
    """Return remainder of a divided by b."""
    if b == 0:
        return "Error: Modulo by zero"
    return a % b

if __name__ == "__main__":
    # Test cases
    print("=== Basic Math Functions ===")
    print(f"Add: 10 + 5 = {add(10, 5)}")
    print(f"Subtract: 10 - 5 = {subtract(10, 5)}")
    print(f"Multiply: 10 * 5 = {multiply(10, 5)}")
    print(f"Divide: 10 / 5 = {divide(10, 5)}")
    print(f"Power: 2 ^ 8 = {power(2, 8)}")
    print(f"Modulo: 17 % 5 = {modulo(17, 5)}")
    print(f"Divide by zero: {divide(10, 0)}")
    print(f"Modulo by zero: {modulo(10, 0)}")
