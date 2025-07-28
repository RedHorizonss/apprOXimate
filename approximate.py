# approximate class for charge balancing chemical formulas
# This class reads a CSV file containing chemical data, parses chemical formulas, and balances charges.
# It supports both fixed and variable oxidation states, and provides detailed debug logging when verbose mode is enabled.
# Translated from the original C++ code into python with now added functionality for interactive use in Jupyter notebooks. see approximate_widgets.ipynb for interactive use.
# Fay Timen and Thomas Ashton 2025 python version 3.11

import pandas as pd
import re
from collections import defaultdict
from decimal import Decimal, getcontext, ROUND_HALF_UP
from fractions import Fraction

class ApprOXimate:
    def __init__(self, csv_file_path = "variable_ox_states_srps.csv", fixed_ox_states_path="fixed_ox_states.csv", verbose=False, precision=6):
        """
        Initialize the ApprOXimate class.
        
        Args:
            csv_file_path (str): Path to the CSV file containing chemical data
            verbose (bool): Enable debug logging
            precision (int): Number of decimal places for precision (default: 6)
            fixed_ox_states_path (str): Path to CSV file with fixed oxidation states (default: "fixed_ox_states.csv")
        """
        self.verbose = verbose
        self.precision = precision
        # Set decimal context for consistent precision
        getcontext().prec = precision + 2  # Extra precision for intermediate calculations
        getcontext().rounding = ROUND_HALF_UP
        self.load_chemical_data(csv_file_path, fixed_ox_states_path)
    
    def precise_round(self, value, decimals=None):
        """
        Round value to specified decimal places using Decimal for precision.
        
        Args:
            value: Value to round (can be float, int, Decimal, or Fraction)
            decimals: Number of decimal places (uses self.precision if None)
        
        Returns:
            float: Precisely rounded value
        """
        if decimals is None:
            decimals = self.precision
        
        if isinstance(value, (int, float)):
            # Convert to Decimal for precise rounding
            decimal_value = Decimal(str(value))
        elif isinstance(value, Fraction):
            # Convert fraction to decimal for rounding
            decimal_value = Decimal(value.numerator) / Decimal(value.denominator)
        elif isinstance(value, Decimal):
            decimal_value = value
        else:
            decimal_value = Decimal(str(value))
        
        # Round using Decimal precision
        rounded_decimal = decimal_value.quantize(Decimal('0.' + '0' * decimals))
        return float(rounded_decimal)
    
    def safe_float_comparison(self, a, b, tolerance=None):
        """
        Compare two floating point numbers with tolerance for precision errors.
        
        Args:
            a, b: Values to compare
            tolerance: Comparison tolerance (uses 10^(-precision) if None)
        
        Returns:
            bool: True if values are equal within tolerance
        """
        if tolerance is None:
            tolerance = 10 ** (-self.precision)
        
        return abs(float(a) - float(b)) <= tolerance
    
    def log(self, message):
        """Print message if verbose mode is enabled"""
        if self.verbose:
            print(f"[DEBUG] {message}")
    
    def load_chemical_data(self, csv_file_path, fixed_ox_states_path):
        self.log(f"Loading chemical data from {csv_file_path}")
        df = pd.read_csv(csv_file_path)
        self.list1 = defaultdict(list)
        self.list2 = {}
        
        for _, row in df.iterrows():
            element = row['Element']
            ox = int(row['Ox'])
            red = int(row['Red'])
            srp = self.precise_round(float(row['SRP (V)']))
            self.list1[element].append((ox, red, srp))
        
        # Sort each element's states by red value (ascending)
        for element in self.list1:
            self.list1[element].sort(key=lambda x: x[1])  # Sort by red value
        
        self.log(f"Loaded {len(self.list1)} elements with variable oxidation states")
        
        # Load fixed oxidation states from CSV file

        self.log(f"Loading fixed oxidation states from {fixed_ox_states_path}")
        fixed_df = pd.read_csv(fixed_ox_states_path)
        
        for _, row in fixed_df.iterrows():
            element = row['Element']
            ox_state = int(row['Ox State'])
            self.list2[element] = ox_state
            self.log(f"Loaded fixed oxidation state: {element} = {ox_state}")
        
        self.log(f"Loaded {len(self.list2)} fixed oxidation states from CSV")
        
        # Remove fixed elements from variable list
        for element in self.list2.keys():
            if element in self.list1:
                del self.list1[element]
                self.log(f"Removed {element} from variable list (now fixed)")
        
        self.log(f"Fixed oxidation states: {self.list2}")
        self.log(f"Variable elements after cleanup: {list(self.list1.keys())}")
    
    def parse_formula(self, formula_input):
        self.log(f"Parsing formula: {formula_input}")
        pattern = r'([A-Z][a-z]*)(\d*\.?\d*/?\d*\.?\d*)'
        matches = re.findall(pattern, formula_input)
        formula = {}
        
        self.log(f"Regex matches: {matches}")
        
        for element, quantity_str in matches:
            if element not in self.list1 and element not in self.list2:
                print(f'Warning: Element {element} not found in database')
                return None
            
            quantity = 1.0
            if quantity_str:
                if '/' in quantity_str:
                    # Use Fraction for exact fractional representation
                    numerator, denominator = quantity_str.split('/')
                    fraction = Fraction(int(numerator), int(denominator))
                    quantity = self.precise_round(fraction)
                    self.log(f"Fractional quantity for {element}: {numerator}/{denominator} = {quantity}")
                else:
                    quantity = self.precise_round(float(quantity_str))
                    self.log(f"Quantity for {element}: {quantity}")
            else:
                self.log(f"No quantity specified for {element}, using 1.0")
            
            # Accumulate quantities with precise arithmetic
            if element in formula:
                formula[element] = self.precise_round(formula[element] + quantity)
            else:
                formula[element] = quantity
        
        self.log(f"Parsed formula: {formula}")
        return formula

    def initialize_element_states(self, formula):
        """Initialize each element to its lowest red value state (where red > 0)"""
        element_states = {}
        
        for element in formula:
            if element in self.list1:
                # Find the first state where red > 0
                initial_state = None
                initial_index = 0
                
                for i, state in enumerate(self.list1[element]):
                    if state[1] > 0:  # red > 0
                        initial_state = state
                        initial_index = i
                        break
                
                if initial_state is None:
                    # If no state has red > 0, use the first state as fallback
                    initial_state = self.list1[element][0]
                    initial_index = 0
                
                element_states[element] = {
                    'current_index': initial_index,
                    'current_state': initial_state,
                    'all_states': self.list1[element]
                }
                self.log(f"Initialized {element} to state: ox={initial_state[0]}, red={initial_state[1]}, srp={initial_state[2]}")
        
        return element_states

    def calculate_current_charge(self, formula, element_states):
        """Calculate the current total charge based on current oxidation states"""
        sum_non_m = Decimal('0')
        sum_m = Decimal('0')
        
        self.log("Calculating charges:")
        for element, quantity in formula.items():
            if element in self.list2:
                red = self.list2[element]
                charge_contribution = Decimal(str(quantity)) * Decimal(str(red))
                sum_non_m += charge_contribution
                self.log(f"  {element} (fixed): {quantity} × {red} = {float(charge_contribution)}")
            else:
                current_state = element_states[element]['current_state']
                charge_contribution = Decimal(str(quantity)) * Decimal(str(current_state[1]))  # Use red value
                sum_m += charge_contribution
                self.log(f"  {element} (variable): {quantity} × {current_state[1]} = {float(charge_contribution)} (ox={current_state[0]}, red={current_state[1]}, srp={current_state[2]})")
        
        total_charge = self.precise_round(sum_non_m + sum_m)
        self.log(f"Total charge: {float(sum_non_m)} + {float(sum_m)} = {total_charge}")
        return total_charge

    def find_lowest_srp_element(self, formula, element_states):
        """Find the element with the lowest SRP among variable elements (excluding fixed elements)"""
        lowest_srp_element = None
        lowest_srp = float('inf')
        
        self.log("Finding element with lowest SRP:")
        for element, quantity in formula.items():
            # Only consider variable elements (not fixed elements in self.list2)
            if element in element_states and element not in self.list2:
                current_state = element_states[element]['current_state']
                self.log(f"  {element}: SRP = {current_state[2]}")
                if current_state[2] < lowest_srp:
                    lowest_srp = current_state[2]
                    lowest_srp_element = {
                        'element': element,
                        'quantity': quantity,
                        'state': current_state
                    }
        
        if lowest_srp_element:
            self.log(f"Element with lowest SRP: {lowest_srp_element['element']} (SRP = {lowest_srp})")
        else:
            self.log("No element with lowest SRP found")
        
        return lowest_srp_element

    def calculate_charge_with_delta(self, lowest_srp_element, total_charge):
        """Calculate what the charge would be if we fully oxidize the lowest SRP element"""
        if not lowest_srp_element:
            return 0, 0
        
        current_state = lowest_srp_element['state']
        delta_ox_calc = (Decimal(str(current_state[0])) - Decimal(str(current_state[1]))) * Decimal(str(lowest_srp_element['quantity']))
        delta_ox_max = self.precise_round(delta_ox_calc)
        
        charge_with_delta = self.precise_round(delta_ox_max + total_charge)
        
        self.log(f"Delta oxidation max: ({current_state[0]} - {current_state[1]}) × {lowest_srp_element['quantity']} = {delta_ox_max}")
        self.log(f"Charge with delta: {delta_ox_max} + {total_charge} = {charge_with_delta}")
        
        return delta_ox_max, charge_with_delta

    def update_oxidation_state(self, lowest_srp_element, element_states):
        """Update the oxidation state of an element to the next available state in stepwise manner"""
        element = lowest_srp_element['element']
        current_index = element_states[element]['current_index']
        all_states = element_states[element]['all_states']
        
        # Move to the next state in the stepwise progression
        if current_index < len(all_states) - 1:
            next_index = current_index + 1
            next_state = all_states[next_index]
            
            element_states[element]['current_index'] = next_index
            element_states[element]['current_state'] = next_state
            
            self.log(f"Updated {element} from index {current_index} to {next_index}")
            self.log(f"New state: ox={next_state[0]}, red={next_state[1]}, srp={next_state[2]}")
        else:
            # If we're at the last state, create a new state where red = ox
            current_state = element_states[element]['current_state']
            new_state = (current_state[0], current_state[0], float('inf'))
            
            # Add the new state to the list and update
            element_states[element]['all_states'].append(new_state)
            element_states[element]['current_index'] = len(element_states[element]['all_states']) - 1
            element_states[element]['current_state'] = new_state
            
            self.log(f"Created new final state for {element}: ox={new_state[0]}, red={new_state[1]}, srp={new_state[2]}")
            
    def charge_balance(self, formula, return_dict=False):
        """
        Main charge balancing function
        
        Args:
            formula: The chemical formula to balance
            return_dict (bool): If True, return results as dictionary; if False, return as string
        
        Returns:
            dict or str: Balanced chemical formula results
        """
        self.log("=== Starting charge balance calculation ===")
        self.log(f"Input formula: {formula}")
        
        # Initialize element states
        element_states = self.initialize_element_states(formula)
        
        # Perform iterative charge balancing
        balance_result = self.perform_charge_balance_iterations(formula, element_states)
        
        # Calculate final charge
        final_charge = self.calculate_final_charge(formula, element_states, balance_result)
        
        # Build result based on return format preference
        if return_dict:
            result = self.build_result_dictionary(formula, element_states, balance_result, final_charge)
        else:
            result = self.build_result_string(formula, element_states, balance_result, final_charge)
        
        return result

    def calculate_final_charge(self, formula, element_states, balance_result):
        """Calculate the final charge after balancing"""
        self.log("\n=== Calculating final charge ===")
        
        element_with_final_lowest_srp = balance_result['element_with_final_lowest_srp']
        final_delta_ox_max = balance_result['final_delta_ox_max']
        final_charge_with_delta = balance_result['final_charge_with_delta']
        element_quantity = balance_result['element_quantity']
        final_ox = balance_result['final_ox']
        final_red = balance_result['final_red']
        
        final_charge = Decimal('0')
        
        # Calculate charge from fixed and other variable elements
        for element, quantity in formula.items():
            if element in self.list2:
                charge_contribution = Decimal(str(quantity)) * Decimal(str(self.list2[element]))
                final_charge += charge_contribution
                self.log(f"Final charge from {element}: {quantity} × {self.list2[element]} = {float(charge_contribution)}")
            elif element != element_with_final_lowest_srp:
                current_state = element_states[element]['current_state']
                charge_contribution = Decimal(str(quantity)) * Decimal(str(current_state[1]))  # Use red value
                final_charge += charge_contribution
                self.log(f"Final charge from {element}: {quantity} × {current_state[1]} = {float(charge_contribution)}")
        
        # Handle the element with final lowest SRP
        if element_with_final_lowest_srp and (self.safe_float_comparison(final_charge_with_delta, 0) or final_charge_with_delta > 0):
            if not self.safe_float_comparison(final_delta_ox_max, 0):  # Only adjust if there's actually a difference in oxidation states
                if final_delta_ox_max != 0:  # Avoid division by zero
                    adjusted_calc = Decimal(str(element_quantity)) * (Decimal(str(final_charge_with_delta)) / Decimal(str(final_delta_ox_max)))
                    adjusted_quant = self.precise_round(adjusted_calc)
                    remaining_quant = self.precise_round(element_quantity - adjusted_quant)
                    
                    self.log(f"Adjusting {element_with_final_lowest_srp}:")
                    self.log(f"  Adjusted quantity: {element_quantity} × ({final_charge_with_delta} / {final_delta_ox_max}) = {adjusted_quant}")
                    self.log(f"  Remaining quantity: {element_quantity} - {adjusted_quant} = {remaining_quant}")
                    
                    if adjusted_quant > 0:
                        charge_contribution = Decimal(str(adjusted_quant)) * Decimal(str(final_red))
                        final_charge += charge_contribution
                        self.log(f"  Charge from adjusted: {adjusted_quant} × {final_red} = {float(charge_contribution)}")
                    if remaining_quant > 0:
                        charge_contribution = Decimal(str(remaining_quant)) * Decimal(str(final_ox))
                        final_charge += charge_contribution
                        self.log(f"  Charge from remaining: {remaining_quant} × {final_ox} = {float(charge_contribution)}")
                else:
                    # If delta is zero, use original quantity
                    charge_contribution = Decimal(str(element_quantity)) * Decimal(str(final_red))
                    final_charge += charge_contribution
                    self.log(f"Delta is zero for {element_with_final_lowest_srp}: {element_quantity} × {final_red} = {float(charge_contribution)}")
            else:
                # If no difference in oxidation states, just use the original quantity
                charge_contribution = Decimal(str(element_quantity)) * Decimal(str(final_red))
                final_charge += charge_contribution
                self.log(f"No oxidation state difference for {element_with_final_lowest_srp}: {element_quantity} × {final_red} = {float(charge_contribution)}")
        
        final_charge = self.precise_round(final_charge)
        self.log(f"Final total charge: {final_charge}")
        
        return final_charge

    def perform_charge_balance_iterations(self, formula, element_states):
        """Perform the iterative charge balancing algorithm"""
        self.log("=== Starting charge balance iterations ===")
        
        total_charge = 0
        previous_charge_with_delta = None  # Use None instead of NaN
        element_with_final_lowest_srp = None
        final_delta_ox_max = 0
        final_charge_with_delta = 0
        element_quantity = 0
        final_ox = 0
        final_red = 0
        
        iteration_count = 0
        max_iterations = 100
        
        while True:
            iteration_count += 1
            self.log(f"\n--- Iteration {iteration_count} ---")
            
            if iteration_count > max_iterations:
                self.log("Maximum iterations reached, breaking")
                break
            
            # Calculate current charge
            total_charge = self.calculate_current_charge(formula, element_states)
            
            # Find element with lowest SRP
            lowest_srp_element = self.find_lowest_srp_element(formula, element_states)
            
            if not lowest_srp_element:
                self.log("No element with lowest SRP found, breaking")
                break
            
            # Calculate charge with delta
            delta_ox_max, charge_with_delta = self.calculate_charge_with_delta(lowest_srp_element, total_charge)
            
            # Use safe comparison for floating point equality
            if previous_charge_with_delta is not None and self.safe_float_comparison(charge_with_delta, previous_charge_with_delta):
                self.log("Charge with delta equals previous, breaking")
                break
            
            previous_charge_with_delta = charge_with_delta
            
            if charge_with_delta < 0:
                self.update_oxidation_state(lowest_srp_element, element_states)
                total_charge = self.precise_round(total_charge + delta_ox_max)
            else:
                self.log("Charge with delta >= 0, finalizing")
                element_with_final_lowest_srp = lowest_srp_element['element']
                final_delta_ox_max = delta_ox_max
                final_charge_with_delta = charge_with_delta
                element_quantity = lowest_srp_element['quantity']
                final_ox = lowest_srp_element['state'][0]
                final_red = lowest_srp_element['state'][1]
                self.log(f"Final element: {element_with_final_lowest_srp}")
                self.log(f"Final delta ox max: {final_delta_ox_max}")
                self.log(f"Final charge with delta: {final_charge_with_delta}")
                break
        
        return {
            'element_with_final_lowest_srp': element_with_final_lowest_srp,
            'final_delta_ox_max': final_delta_ox_max,
            'final_charge_with_delta': final_charge_with_delta,
            'element_quantity': element_quantity,
            'final_ox': final_ox,
            'final_red': final_red
        }
    
    def build_result_dictionary(self, formula, element_states, balance_result, final_charge):
        """
        Build the final result as a dictionary
        
        Returns:
            dict: Dictionary containing element information and final charge balance
        """
        self.log("\n=== Building result dictionary ===")
        
        element_with_final_lowest_srp = balance_result['element_with_final_lowest_srp']
        final_delta_ox_max = balance_result['final_delta_ox_max']
        final_charge_with_delta = balance_result['final_charge_with_delta']
        element_quantity = balance_result['element_quantity']
        final_ox = balance_result['final_ox']
        final_red = balance_result['final_red']
        
        result_dict = {
            'elements': {},
            'final_charge_balance': self.precise_round(final_charge, 3)
        }
        
        # Add fixed elements first
        for element in sorted(formula.keys()):
            if element in self.list2:
                quantity = self.precise_round(formula[element], 3)
                result_dict['elements'][element] = {
                    'oxidation_state': self.list2[element],
                    'quantity': quantity,
                    'type': 'fixed'
                }
                self.log(f"Added fixed element to dict: {element} = {result_dict['elements'][element]}")
        
        # Add variable elements
        for element in sorted(formula.keys()):
            if element in element_states and element != element_with_final_lowest_srp:
                current_state = element_states[element]['current_state']
                quantity = self.precise_round(formula[element], 3)
                result_dict['elements'][element] = {
                    'oxidation_state': current_state[1],  # Using red value
                    'quantity': quantity,
                    'type': 'variable',
                    'srp': current_state[2]
                }
                self.log(f"Added variable element to dict: {element} = {result_dict['elements'][element]}")
        
        # Handle element with lowest SRP (may have split quantities)
        if element_with_final_lowest_srp and (self.safe_float_comparison(final_charge_with_delta, 0) or final_charge_with_delta > 0):
            if not self.safe_float_comparison(final_delta_ox_max, 0):
                if final_delta_ox_max != 0:  # Avoid division by zero
                    adjusted_calc = Decimal(str(element_quantity)) * (Decimal(str(final_charge_with_delta)) / Decimal(str(final_delta_ox_max)))
                    adjusted_quant = self.precise_round(adjusted_calc, 3)
                    remaining_quant = self.precise_round(element_quantity - adjusted_quant, 3)
                    
                    # Create a list to store multiple oxidation states for this element
                    element_states_list = []
                    
                    if adjusted_quant > 0:
                        element_states_list.append({
                            'oxidation_state': final_red,
                            'quantity': adjusted_quant
                        })
                        self.log(f"Added adjusted lowest SRP element to dict: {element_with_final_lowest_srp} (red state)")
                    
                    if remaining_quant > 0:
                        element_states_list.append({
                            'oxidation_state': final_ox,
                            'quantity': remaining_quant
                        })
                        self.log(f"Added remaining lowest SRP element to dict: {element_with_final_lowest_srp} (ox state)")
                    
                    result_dict['elements'][element_with_final_lowest_srp] = {
                        'states': element_states_list,
                        'total_quantity': self.precise_round(element_quantity, 3),
                        'type': 'variable_split',
                        'srp': self.precise_round([state[2] for state in element_states[element_with_final_lowest_srp]['all_states'] if state[1] == final_red][0], 3)
                    }
                else:
                    # If delta is zero, use original quantity
                    quantity = self.precise_round(element_quantity, 3)
                    result_dict['elements'][element_with_final_lowest_srp] = {
                        'oxidation_state': final_red,
                        'quantity': quantity,
                        'type': 'variable',
                        'srp': self.precise_round([state[2] for state in element_states[element_with_final_lowest_srp]['all_states'] if state[1] == final_red][0], 3)
                    }
                    self.log(f"Added lowest SRP element to dict (delta zero): {element_with_final_lowest_srp}")
            else:
                # If no difference in oxidation states, just show the original quantity
                quantity = self.precise_round(element_quantity, 3)
                result_dict['elements'][element_with_final_lowest_srp] = {
                    'oxidation_state': final_red,
                    'quantity': quantity,
                    'type': 'variable',
                    'srp': self.precise_round([state[2] for state in element_states[element_with_final_lowest_srp]['all_states'] if state[1] == final_red][0], 3)
                }
                self.log(f"Added lowest SRP element to dict (no difference): {element_with_final_lowest_srp}")
        
        self.log(f"Final result dictionary: {result_dict}")
        return result_dict
    
    def charge_balance_dict(self, formula):
        """
        Convenience method that returns results as a dictionary
        
        Args:
            formula: The chemical formula to balance
        
        Returns:
            dict: Balanced chemical formula results in dictionary format
        """
        return self.charge_balance(formula, return_dict=True)
    
    def build_result_string(self, formula, element_states, balance_result, final_charge):
        """Build the final result string"""
        self.log("\n=== Building result string ===")
        
        element_with_final_lowest_srp = balance_result['element_with_final_lowest_srp']
        final_delta_ox_max = balance_result['final_delta_ox_max']
        final_charge_with_delta = balance_result['final_charge_with_delta']
        element_quantity = balance_result['element_quantity']
        final_ox = balance_result['final_ox']
        final_red = balance_result['final_red']
        
        result_parts = []
        
        # Add fixed elements first
        for element in sorted(formula.keys()):
            if element in self.list2:
                quantity = self.precise_round(formula[element], 3)  # Round to 3 decimals for display
                part = f"{element}:{self.list2[element]}:{quantity}"
                result_parts.append(part)
                self.log(f"Added fixed element: {part}")
        
        # Add variable elements
        for element in sorted(formula.keys()):
            if element in element_states and element != element_with_final_lowest_srp:
                current_state = element_states[element]['current_state']
                quantity = self.precise_round(formula[element], 3)  # Round to 3 decimals for display
                part = f"{element}:{current_state[1]}:{quantity}"
                result_parts.append(part)
                self.log(f"Added variable element: {part}")
        
        # Handle element with lowest SRP
        if element_with_final_lowest_srp and (self.safe_float_comparison(final_charge_with_delta, 0) or final_charge_with_delta > 0):
            if not self.safe_float_comparison(final_delta_ox_max, 0):
                if final_delta_ox_max != 0:  # Avoid division by zero
                    adjusted_calc = Decimal(str(element_quantity)) * (Decimal(str(final_charge_with_delta)) / Decimal(str(final_delta_ox_max)))
                    adjusted_quant = self.precise_round(adjusted_calc, 3)
                    remaining_quant = self.precise_round(element_quantity - adjusted_quant, 3)
                    
                    if adjusted_quant > 0:
                        part = f"{element_with_final_lowest_srp}:{final_red}:{adjusted_quant}"
                        result_parts.append(part)
                        self.log(f"Added adjusted lowest SRP element: {part}")
                    if remaining_quant > 0:
                        part = f"{element_with_final_lowest_srp}:{final_ox}:{remaining_quant}"
                        result_parts.append(part)
                        self.log(f"Added remaining lowest SRP element: {part}")
                else:
                    # If delta is zero, use original quantity
                    quantity = self.precise_round(element_quantity, 3)
                    part = f"{element_with_final_lowest_srp}:{final_red}:{quantity}"
                    result_parts.append(part)
                    self.log(f"Added lowest SRP element (delta zero): {part}")
            else:
                # If no difference in oxidation states, just show the original quantity
                quantity = self.precise_round(element_quantity, 3)
                part = f"{element_with_final_lowest_srp}:{final_red}:{quantity}"
                result_parts.append(part)
                self.log(f"Added lowest SRP element (no difference): {part}")
        
        final_charge_rounded = self.precise_round(final_charge, 3)
        result_parts.append(f"FinalChargeBalance:{final_charge_rounded}")
        self.log(f"Added final charge balance: FinalChargeBalance:{final_charge_rounded}")
        
        result = ";".join(result_parts)
        self.log(f"Final result: {result}")
        
        return result