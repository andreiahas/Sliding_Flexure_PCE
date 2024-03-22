import os
import pandas as pd

class InputChecker:
    def check_input(T1, muf, mratio, Tb, RQI, IM):
        """
        Check if the input parameters are within specified ranges.

        Args:
        - T1: superstructure first-mode natural period T1.
        - muf: base-isolation friction coefficient
        - mratio: mass ratio m1*/mtot
        - Tb: Base isolation period Tb
        - RQI: response quantity of interest
        - IM: intensity measure

        Return:
        - If any input parameter is outside the specified range, returns an error message.
        Otherwise, returns None.
        """
        if not (0.1 <= T1 <= 1):
            return "Error: T1 must be in the range [0.1, 1]"
    
        if not (0.03 <= muf <= 0.18):
            return "Error: muf must be in the range [0.03, 0.18]"

        if not (0.3 <= mratio <= 0.9):
            return "Error: mratio must be in the range [0.3, 0.9]"

        if not (3 <= Tb <= 6):
            return "Error: Tb must be in the range [3, 6]"
    
        if RQI not in ['D1', 'u0']:
            return "Error: RQI must be 'D1' or 'u0'"

        if IM not in ['SC', 'GM']:
            return "Error: IM must be 'SC' or 'GM'"

        return None  # No errors found

class DataLoader:
    def load_DF_from_CSV_sliding(IM, RQI, var):
        """
        Load DataFrames based on the input parameters.

        Args:
        - IM: Value of IM.
        - RQI: Value of RQI.
        - var: Value of var.

        Return:
        - DataFrame loaded from the CSV file.
        """
        # Obtain the name of the current directory
        current_directory = os.getcwd()

        # Determine load_name based on RQI value
        if RQI == 'D1':
            load_name = 'D1_Sd'
        elif RQI == 'u0':
            load_name = 'u0_PGV'
        else:
            raise ValueError("Invalid value for RQI. Must be 'D1' or 'u0'.")

        # Construct the file name
        file_name_basis = f"CSVfiles\PCE_{var}_{load_name}{IM}_BASIS.csv"
        file_name_coefs = f"CSVfiles\PCE_{var}_{load_name}{IM}_COEFS.csv"

        # Load the DataFrame from the CSV file
        file_path_basis = os.path.join(current_directory, file_name_basis)
        file_path_coefs = os.path.join(current_directory, file_name_coefs)
    
        try:
            df_basis = pd.read_csv(file_path_basis, header=None)
            #print("Basis CSV file loaded successfully.")
        except FileNotFoundError:
            print("Basis file not found. Please make sure the file exists in the current directory.")
        except Exception as e:
            print("An error occurred:", e)
    
        try:
            df_coefs = pd.read_csv(file_path_coefs, header=None)
            #print("Coefs CSV file loaded successfully.")
        except FileNotFoundError:
            print("Coefs file not found. Please make sure the file exists in the current directory.")
        except Exception as e:
            print("An error occurred:", e)
        
        return df_basis, df_coefs

class LegendrePolynomials:
    def map_to_interval_m1to1(Z, a, b):
        """
        Map a value Z within the interval [a, b] to the interval [-1, 1].

        Args:
        - Z: The value to be mapped.
        - a: Lower bound of the interval [a, b].
        - b: Upper bound of the interval [a, b].

        Return:
        - The mapped value within the interval [-1, 1].
        """
        # Map the proportional value to the [-1, 1] interval
        mapped_value = 2 * (Z - a) / (b - a) - 1

        return mapped_value

    def normalized_legendre_polynomials(x_values,keys):
        """
        Calculate normalized Legendre polynomials of degrees 1 to 7 for multiple input values.

        Args:
        - x_values: A list or tuple of input values at which to evaluate the normalized Legendre polynomials.
        - keys: the name of the variables that will be the keys in the dictionary

        Return:
        - A dictionary where keys are input value position and values are dictionaries containing
        normalized Legendre polynomials of degrees 1 to 7 evaluated at the corresponding input value.
        """
        if len(x_values) != len(keys):
            raise ValueError("Error in normalized_legendre_polynomials: x_values and keys must have the same number of elements")
    
        normalized_legendre_values_all = {}

        for i in range(0,len(x_values)):
            x = x_values[i]
            normalized_legendre_values = {}
        
            # Degree 0 normalized Legendre polynomial
            normalized_legendre_values[0] = 1

            # Degree 1 normalized Legendre polynomial
            normalized_legendre_values[1] = x * (3 ** 0.5)

            # Degree 2 normalized Legendre polynomial
            normalized_legendre_values[2] = 0.5 * (3 * x**2 - 1) * (5 ** 0.5)

            # Degree 3 normalized Legendre polynomial
            normalized_legendre_values[3] = 0.5 * (5 * x**3 - 3 * x) * (7 ** 0.5)

            # Degree 4 normalized Legendre polynomial
            normalized_legendre_values[4] = (1/8) * (35 * x**4 - 30 * x**2 + 3) * (9 ** 0.5)

            # Degree 5 normalized Legendre polynomial
            normalized_legendre_values[5] = (1/8) * (63 * x**5 - 70 * x**3 + 15 * x) * (11 ** 0.5)

            # Degree 6 normalized Legendre polynomial
            normalized_legendre_values[6] = (1/16) * (231 * x**6 - 315 * x**4 + 105 * x**2 - 5) * (13 ** 0.5)

            # Degree 7 normalized Legendre polynomial
            normalized_legendre_values[7] = (1/16) * (429 * x**7 - 693 * x**5 + 315 * x**3 - 35 * x) * (15 ** 0.5)

            normalized_legendre_values_all[keys[i]] = normalized_legendre_values

        return normalized_legendre_values_all

class PCECalculator:
    def calculate_PCE(Ldict, df_basis, df_coefs):
        """
        Calculate PCE based on the input data.

        Args:
        - Ldict: Nested dictionary containing legendre polynomials
        - df_basis: DataFrame containing basis values.
        - df_coefs: DataFrame containing coefficients.

        Return:
        - YPCE: Calculated value of PCE.
        """
        # Initialize PCE as a numeric value
        YPCE = 0

        # Iterate over the values of df_coefs with index i
        for i in range(len(df_coefs)):
            # Get coefficients from df_coefs at index i
            coef_i = df_coefs.iloc[i, 0]

            # Get basis values from df_basis at index i
            basis_values_i = df_basis.iloc[i]

            # Calculate product and add to the coefficient
            YPCE += coef_i * Ldict['T1'][basis_values_i[0]] * Ldict['mratio'][basis_values_i[1]] * Ldict['muf'][basis_values_i[2]] * Ldict['Tb'][basis_values_i[3]]

        return YPCE

    def run_PCE(T1,muf,mratio,Tb,RQI,IM):
        """
        Calculate PCE metamodel results for coefficients c1, c2, and beta.
    
        Args:
        - T1: superstructure first-mode natural period T1.
        - muf: base-isolation friction coefficient
        - mratio: mass ratio m1*/mtot
        - Tb: base isolation period Tb
        - RQI: response quantity of interest
        - IM: intensity measure

        Return:
        - c1, c2, beta: Calculated coefficients c1, c2, and beta with the PCE metamodel.
        """
        ## PCE framework:
        # Check if inputs are in the possible range:
        if InputChecker.check_input(T1, muf, mratio, Tb, RQI, IM) is not None:
            print(InputChecker.check_input(T1, muf, mratio, Tb, RQI, IM))
            return None

        # Load data frames with basis and coefficients of PCE:
        df_basis_C1, df_coefs_C1 = DataLoader.load_DF_from_CSV_sliding(IM, RQI, 'C1')
        df_basis_C2, df_coefs_C2 = DataLoader.load_DF_from_CSV_sliding(IM, RQI, 'C2')
        df_basis_beta, df_coefs_beta = DataLoader.load_DF_from_CSV_sliding(IM, RQI, 'beta')

        # Convert the input into the [-1,1] interval:
        values_and_ranges = [
            (T1, 0.1, 1),
            (mratio, 0.3, 0.9),
            (muf, 0.03, 0.18),
            (Tb, 3, 6)
        ]
        keys = ['T1', 'mratio', 'muf', 'Tb']

        # Create the array X from transformation of vector Z using list comprehension:
        X = [LegendrePolynomials.map_to_interval_m1to1(value, low, high) for value, low, high in values_and_ranges]

        # This function returns a dictonary with Legendre polynomials (degree,value), inside of a dictionary linked to the keys
        legendre_dict_X = LegendrePolynomials.normalized_legendre_polynomials(X,keys)

        c1 = PCECalculator.calculate_PCE(legendre_dict_X, df_basis_C1, df_coefs_C1)
        c2 = PCECalculator.calculate_PCE(legendre_dict_X, df_basis_C2, df_coefs_C2)
        beta = PCECalculator.calculate_PCE(legendre_dict_X, df_basis_beta, df_coefs_beta)
    
        return c1, c2, beta