import tkinter as tk
from tkinter import simpledialog, messagebox
import pandas as pd
import numpy as np

# Constants
R = 8.314  # Gas constant in J/(mol·K)

# File paths
BINARY_ENTHALPY_TAKEUCHI = "binary_enthalpies_takeuchi.xlsx"
BINARY_ENTHALPY_TROPAREVSKY = "binary_enthalpies_troparevsky.xlsx"
ATOMIC_RADII_FILE = "atomic_radii.xlsx"

# Periodic Table Layout
PERIODIC_TABLE = [
    ["H", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "He"],
    ["Li", "Be", "", "", "", "", "", "", "", "", "", "", "B", "C", "N", "O", "F", "Ne"],
    ["Na", "Mg", "", "", "", "", "", "", "", "", "", "", "Al", "Si", "P", "S", "Cl", "Ar"],
    ["K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr"],
    ["Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe"],
    ["Cs", "Ba", "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn"],
    ["Fr", "Ra", "Ac", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"],
    ["", "", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"],
    ["", "", "Th", "Pa", "U", "Np", "Pu"]
]

# Load atomic radii data
def load_atomic_radii(file_path):
    df = pd.read_excel(file_path)
    return dict(zip(df["Symbol"], df["Atomic Radius (pm)"]))

atomic_radii_data = load_atomic_radii(ATOMIC_RADII_FILE)

# Load binary enthalpy database
def read_binary_enthalpy_database(file_path):
    return pd.read_excel(file_path, index_col=0)

# Check if elements exist in the binary enthalpy database
def check_database_compatibility(elements, database):
    return all(element in database.index and element in database.columns for element in elements)

# Calculate mixing enthalpy
def calculate_mixing_enthalpy(elements, atomic_percents, binary_database):
    n = len(elements)
    mole_fractions = np.array(atomic_percents) / np.sum(atomic_percents)
    mixing_enthalpy = sum(4 * mole_fractions[i] * mole_fractions[j] * binary_database.loc[elements[i], elements[j]]
                          for i in range(n) for j in range(i + 1, n))
    return mixing_enthalpy

# Calculate mixing entropy
def calculate_mixing_entropy(atomic_percents):
    mole_fractions = np.array(atomic_percents) / np.sum(atomic_percents)
    return -R * np.sum(mole_fractions * np.log(mole_fractions))

# Calculate free energy
def calculate_mixing_free_energy(mixing_enthalpy, mixing_entropy, temperature):
    return mixing_enthalpy - (temperature * (mixing_entropy / 1000))

# Calculate delta parameter
def calculate_delta_parameter(elements, atomic_percents):
    mole_fractions = np.array(atomic_percents) / np.sum(atomic_percents)
    radii = np.array([atomic_radii_data[el] for el in elements])

    avg_radius = np.sum(mole_fractions * radii)
    delta = np.sqrt(np.sum(mole_fractions * ((radii - avg_radius) / avg_radius) ** 2))
    return delta

# GUI Application
class PeriodicTableApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Multielement Alloy Free Energy Calculator")

        self.selected_elements = []
        self.atomic_percents = []
        self.reaction_temperature = 293.15  # Default temperature

        self.create_periodic_table()
        self.create_control_buttons()

    def create_periodic_table(self):
        self.buttons = {}
        frame = tk.Frame(self.root)
        frame.pack()

        for r, row in enumerate(PERIODIC_TABLE):
            for c, element in enumerate(row):
                if element:
                    btn = tk.Button(frame, text=element, width=5, height=2,
                                    command=lambda e=element: self.select_element(e))
                    btn.grid(row=r, column=c)
                    self.buttons[element] = btn

    def create_control_buttons(self):
        frame = tk.Frame(self.root)
        frame.pack(pady=10)

        tk.Label(frame, text="Reaction Temperature (K):").pack(side=tk.LEFT, padx=5)
        self.temp_entry = tk.Entry(frame, width=10)
        self.temp_entry.insert(0, "293.15")
        self.temp_entry.pack(side=tk.LEFT, padx=5)

        tk.Button(frame, text="Calculate", command=self.calculate).pack(side=tk.LEFT, padx=5)
        tk.Button(frame, text="Reset", command=self.reset).pack(side=tk.LEFT, padx=5)

    def select_element(self, element):
        if element in self.selected_elements:
            messagebox.showwarning("Warning", f"{element} is already selected.")
            return

        remaining_percentage = 100 - sum(self.atomic_percents)
        percent = simpledialog.askfloat("Input", f"Enter atomic percent for {element} (Remaining: {remaining_percentage:.2f}%):")

        if percent is None or percent <= 0 or sum(self.atomic_percents) + percent > 100:
            messagebox.showerror("Error", "Invalid percentage input.")
            return

        self.selected_elements.append(element)
        self.atomic_percents.append(percent)
        self.buttons[element].config(bg="lightblue", state=tk.DISABLED)

    def calculate(self):
        if abs(sum(self.atomic_percents) - 100) > 1e-6:
            messagebox.showerror("Error", "Atomic percentages must sum to 100%.")
            return

        try:
            self.reaction_temperature = float(self.temp_entry.get())
        except ValueError:
            messagebox.showerror("Error", "Invalid temperature input.")
            return

        db_takeuchi = read_binary_enthalpy_database(BINARY_ENTHALPY_TAKEUCHI)
        db_troparevsky = read_binary_enthalpy_database(BINARY_ENTHALPY_TROPAREVSKY)

        takeuchi_compatible = check_database_compatibility(self.selected_elements, db_takeuchi)
        troparevsky_compatible = check_database_compatibility(self.selected_elements, db_troparevsky)

        mixing_entropy = calculate_mixing_entropy(self.atomic_percents)
        delta = calculate_delta_parameter(self.selected_elements, self.atomic_percents)

        results = f"Delta Parameter: {delta:.4f}\nMixing Entropy: {mixing_entropy:.4f} J/K·mol\n"

        for name, db, compat in [("Takeuchi", db_takeuchi, takeuchi_compatible),
                                 ("Troparevsky", db_troparevsky, troparevsky_compatible)]:
            if compat:
                enthalpy = calculate_mixing_enthalpy(self.selected_elements, self.atomic_percents, db)
                results += (f"\n{name} Mixing Enthalpy: {enthalpy:.4f} kJ/mol\n"
                            f"Free Energy (293.15 K): {calculate_mixing_free_energy(enthalpy, mixing_entropy, 293.15):.4f} kJ/mol\n"
                            f"Free Energy ({self.reaction_temperature:.2f} K): {calculate_mixing_free_energy(enthalpy, mixing_entropy, self.reaction_temperature):.4f} kJ/mol\n")

        messagebox.showinfo("Results", results)

    def reset(self):
        self.selected_elements, self.atomic_percents = [], []
        for btn in self.buttons.values():
            btn.config(bg="SystemButtonFace", state=tk.NORMAL)

# Run GUI
if __name__ == "__main__":
    root = tk.Tk()
    app = PeriodicTableApp(root)
    root.mainloop()
