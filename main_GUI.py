import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np
from backend import get_ortho_expressions
from sympy import integrate, Symbol, factor, exp, oo, sin, pi, cos, Abs, sympify, lambdify, sqrt, powsimp
import scipy.integrate as spi
import os
import pickle
import threading

x = Symbol('x', real=True)

poly_base_params = {"Hermite":(-oo, oo, exp(-(x**2))),
                    "Laguerre":(0, oo, exp(-x)),
                    "Legendre":(-1, 1, 1),
                    "Chebyshev":(-1, 1, sqrt(1 - x**2))
                    }

default_functions = {"sin(x)" : sin(x),
                     "cos(x)": cos(x),
                     "exp(x)" : exp(x),
                     "Abs(x)": Abs(x)}

class OrthogonalPolyGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Orthogonal Polynomial Visualizer")
        self.geometry("1000x800")
        self.minsize(1000, 800)

        # Configure grid to allow resizing
        self.columnconfigure(0, weight=1) # Control panel
        self.columnconfigure(1, weight=4) # Plot area
        self.rowconfigure(0, weight=1)

        # ==========================================
        # MULTI-THREADED PERSISTENT CACHE
        # ==========================================
        self.cache_file = "orthogonal_polys_cache.pkl"
        self.cache_lock = threading.Lock()
        self.basis_cache = self.load_cache()
        
        self.max_calc_degree = 10
        
        self.cached_func_str = ""
        self.cached_coefs = []
        self.last_plotted_degree = -1

        self.create_widgets()
        self.update_ui_state()
        
        # Start the background pre-calculation thread quietly!
        # "daemon=True" means this thread will automatically die when you close the GUI.
        threading.Thread(target=self.background_precalc, daemon=True).start()
        # ==========================================
        # Intercept the window closing event to save the cache
        self.protocol("WM_DELETE_WINDOW", self.on_closing)

    def on_closing(self):
        """Saves the cache to the hard drive right before the application closes."""
        print("Saving polynomial cache to disk... please wait...")
        self.save_cache()
        print("Save complete. Goodbye!")
        self.destroy()

    def load_cache(self):
        """Loads the saved polynomials from the hard drive if the file exists."""
        if os.path.exists(self.cache_file):
            try:
                with open(self.cache_file, "rb") as f:
                    print("Loaded saved polynomial cache from disk!")
                    return pickle.load(f)
            except Exception as e:
                print(f"Failed to load cache: {e}")
        return {}

    def save_cache(self):
        """Safely writes the current dictionary to the hard drive."""
        with self.cache_lock:
            try:
                with open(self.cache_file, "wb") as f:
                    pickle.dump(self.basis_cache, f)
            except Exception as e:
                print(f"Failed to save cache: {e}")

    def background_precalc(self):
        """Quietly calculates the 4 standard bases in the background."""
        for basis_name, (a, b, sigma) in poly_base_params.items():
            # If it's not in the cache, or if we don't have enough polynomials...
            if basis_name not in self.basis_cache:
                print(f"[Background Task] Crunching Gram-Schmidt for {basis_name}...")

                polys = get_ortho_expressions(self.max_calc_degree + 1, a, b, sigma)
                
                with self.cache_lock:
                    self.basis_cache[basis_name] = polys
                    
                self.save_cache()
                print(f"[Background Task] {basis_name} complete and saved!")

            elif len(self.basis_cache[basis_name]) < self.max_calc_degree + 1:
                print(f"[Background Task] Crunching More Gram-Schmidt for {basis_name}...")
                
                polys = get_ortho_expressions(self.max_calc_degree + 1, a, b, sigma, start=self.basis_cache[basis_name])
                
                with self.cache_lock:
                    self.basis_cache[basis_name] = polys
                    
                self.save_cache()
                print(f"[Background Task] {basis_name} complete and saved!")
                
        # print("[Background Task] All standard bases are ready to go!")    

    def create_widgets(self):
        # --- Control Panel (Left Side) ---
        control_frame = ttk.LabelFrame(self, text="Settings", padding="15")
        control_frame.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)

        # 1. Operation Mode
        ttk.Label(control_frame, text="Operation Mode:").pack(anchor="w", pady=(0, 5))
        self.mode_var = tk.StringVar(value="approx")
        self.mode_var.trace_add("write", self.update_ui_state) # Trigger UI update on change
        
        ttk.Radiobutton(control_frame, text="Function Approximation", variable=self.mode_var, value="approx").pack(anchor="w")
        ttk.Radiobutton(control_frame, text="Display First 'n' Polynomials", variable=self.mode_var, value="basis").pack(anchor="w", pady=(0, 15))

        # 2. Polynomial Basis Section
        self.basis_frame = ttk.LabelFrame(control_frame, text="Polynomial Basis", padding="10")
        self.basis_frame.pack(fill="x", pady=(0, 15))

        self.basis_input_type_var = tk.StringVar(value="preset")
        self.basis_input_type_var.trace_add("write", self.update_ui_state)

        # 2a. Preset Basis
        ttk.Radiobutton(self.basis_frame, text="Preset Basis", variable=self.basis_input_type_var, value="preset").pack(anchor="w")
        self.basis_var = tk.StringVar(value="Legendre")
        self.basis_dropdown = ttk.Combobox(self.basis_frame, textvariable=self.basis_var, state="readonly")
        self.basis_dropdown['values'] = tuple(poly_base_params.keys())
        self.basis_dropdown.pack(fill="x", padx=20, pady=(0, 10))
        self.basis_var.trace_add("write", self.on_preset_basis_change)

        # 2b. Custom Basis
        ttk.Radiobutton(self.basis_frame, text="Custom Basis", variable=self.basis_input_type_var, value="custom").pack(anchor="w")
        
        custom_basis_container = ttk.Frame(self.basis_frame)
        custom_basis_container.pack(fill="x", padx=20, pady=(0, 5))
        
        ttk.Label(custom_basis_container, text="a:").grid(row=0, column=0, sticky="w")
        self.custom_a_var = tk.StringVar(value="-1")
        self.custom_a_entry = ttk.Entry(custom_basis_container, textvariable=self.custom_a_var, width=7)
        self.custom_a_entry.grid(row=0, column=1, padx=(0, 10))

        ttk.Label(custom_basis_container, text="b:").grid(row=0, column=2, sticky="w")
        self.custom_b_var = tk.StringVar(value="1")
        self.custom_b_entry = ttk.Entry(custom_basis_container, textvariable=self.custom_b_var, width=7)
        self.custom_b_entry.grid(row=0, column=3)

        ttk.Label(custom_basis_container, text="sigma(x):").grid(row=1, column=0, columnspan=4, sticky="w", pady=(5,0))
        self.custom_sigma_var = tk.StringVar(value="1")
        self.custom_sigma_entry = ttk.Entry(custom_basis_container, textvariable=self.custom_sigma_var)
        self.custom_sigma_entry.grid(row=2, column=0, columnspan=4, sticky="ew")
        custom_basis_container.columnconfigure(3, weight=1)

        # 3. Function Input Section
        self.func_frame = ttk.LabelFrame(control_frame, text="Target Function", padding="10")
        self.func_frame.pack(fill="x", pady=(0, 15))

        self.input_type_var = tk.StringVar(value="preset")
        self.input_type_var.trace_add("write", self.update_ui_state)

        # 3a. Preset Functions
        ttk.Radiobutton(self.func_frame, text="Preset Function", variable=self.input_type_var, value="preset").pack(anchor="w")
        self.preset_var = tk.StringVar(value="sin(x)")
        self.preset_dropdown = ttk.Combobox(self.func_frame, textvariable=self.preset_var, state="readonly")
        self.preset_dropdown['values'] = list(default_functions.values())
        self.preset_dropdown.pack(fill="x", padx=20, pady=(0, 10))
        self.preset_var.trace_add("write", lambda *args: self.on_plot_click())

        # 3b. Custom Function
        ttk.Radiobutton(self.func_frame, text="Custom Function f(x)", variable=self.input_type_var, value="custom").pack(anchor="w")
        self.custom_func_var = tk.StringVar(value="x**2 + 2*x")
        self.custom_entry = ttk.Entry(self.func_frame, textvariable=self.custom_func_var)
        self.custom_entry.pack(fill="x", padx=20, pady=(0, 5))

        # 4. Slider for Approximation Degree / Number of Polynomials
        self.degree_label_text = tk.StringVar(value="Degree / Max 'n':")
        ttk.Label(control_frame, textvariable=self.degree_label_text).pack(anchor="w", pady=(5, 5))
        
        self.degree_var = tk.DoubleVar(value=5)
        slider_frame = ttk.Frame(control_frame)
        slider_frame.pack(fill="x", pady=(0, 20))
        
        self.degree_slider = ttk.Scale(
            slider_frame, from_=1, to=self.max_calc_degree, orient="horizontal",
            variable=self.degree_var, command=self.update_slider_label
        )
        self.degree_slider.pack(side="left", fill="x", expand=True)
        
        self.degree_val_label = ttk.Label(slider_frame, text="5")
        self.degree_val_label.pack(side="right", padx=(10, 0))

        # 5. Axis Limits
        axis_frame = ttk.LabelFrame(control_frame, text="Axis Limits (Leave blank for Auto Y)", padding="10")
        axis_frame.pack(fill="x", pady=(0, 15))

        # X Limits (Side by side)
        x_lim_frame = ttk.Frame(axis_frame)
        x_lim_frame.pack(fill="x", pady=2)
        ttk.Label(x_lim_frame, text="X Min:").pack(side="left")
        self.x_min_var = tk.StringVar(value="-1.0")
        ttk.Entry(x_lim_frame, textvariable=self.x_min_var, width=8).pack(side="left", padx=(0, 10))
        
        ttk.Label(x_lim_frame, text="X Max:").pack(side="left")
        self.x_max_var = tk.StringVar(value="1.0")
        ttk.Entry(x_lim_frame, textvariable=self.x_max_var, width=8).pack(side="left")

        # Y Limits (Side by side)
        y_lim_frame = ttk.Frame(axis_frame)
        y_lim_frame.pack(fill="x", pady=2)
        ttk.Label(y_lim_frame, text="Y Min:").pack(side="left")
        self.y_min_var = tk.StringVar(value="")
        ttk.Entry(y_lim_frame, textvariable=self.y_min_var, width=8).pack(side="left", padx=(0, 10))
        
        ttk.Label(y_lim_frame, text="Y Max:").pack(side="left")
        self.y_max_var = tk.StringVar(value="")
        ttk.Entry(y_lim_frame, textvariable=self.y_max_var, width=8).pack(side="left")

        # 6. Update Button
        self.update_btn = ttk.Button(
            control_frame, 
            text="Plot Approximation", 
            command=self.on_plot_click
        )
        self.update_btn.pack(fill="x", pady=(10, 0))

        # --- Plot Area (Right Side) ---
        plot_frame = ttk.Frame(self, padding="10")
        plot_frame.grid(row=0, column=1, sticky="nsew")

        self.figure = Figure(figsize=(6, 5), dpi=100)
        self.ax = self.figure.add_subplot(111)
        
        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        self.on_plot_click()

    def sync_preset_to_custom(self):
        """Fills the custom basis text boxes with the currently selected preset values."""
        basis_name = self.basis_var.get()
        if basis_name in poly_base_params:
            a, b, sigma = poly_base_params[basis_name]
            
            # str() beautifully converts SymPy objects like 'oo' into plain text strings
            self.custom_a_var.set(str(a))
            self.custom_b_var.set(str(b))
            self.custom_sigma_var.set(str(sigma))

    def on_preset_basis_change(self, *args):
        """Handles when the user picks a new preset basis from the dropdown."""
        # Only overwrite the custom boxes if we are actually in preset mode
        if self.basis_input_type_var.get() == "preset":
            self.sync_preset_to_custom()
        self.on_plot_click()

    def update_ui_state(self, *args):
        """Enables or disables widgets dynamically based on user selections."""
        mode = self.mode_var.get()
        input_type = self.input_type_var.get()

        if mode == "basis":
            # Disable the entire function frame if we are just showing basis polynomials
            for child in self.func_frame.winfo_children():
                child.configure(state="disabled")
            self.degree_label_text.set("Show first 'n' polynomials:")
        else:
            # Enable the radio buttons inside the function frame
            for child in self.func_frame.winfo_children():
                # Re-enable radio buttons
                if isinstance(child, ttk.Radiobutton):
                    child.configure(state="normal")
            self.degree_label_text.set("Approximation Degree (n):")

            # Toggle preset vs custom inputs based on radio selection
            if input_type == "preset":
                self.preset_dropdown.configure(state="readonly")
                self.custom_entry.configure(state="disabled")
            else:
                self.preset_dropdown.configure(state="disabled")
                self.custom_entry.configure(state="normal")

        # Toggle preset vs custom basis inputs
        basis_input_type = self.basis_input_type_var.get()
        if basis_input_type == "preset":
            self.basis_dropdown.configure(state="readonly")
            self.custom_a_entry.configure(state="disabled")
            self.custom_b_entry.configure(state="disabled")
            self.custom_sigma_entry.configure(state="disabled")
            self.sync_preset_to_custom()
        else:
            self.basis_dropdown.configure(state="disabled")
            self.custom_a_entry.configure(state="normal")
            self.custom_b_entry.configure(state="normal")
            self.custom_sigma_entry.configure(state="normal")

        self.on_plot_click()

    def update_slider_label(self, event):
        # The slider internally uses a DoubleVar (e.g., 4.123, 4.124...)
        # We snap it to an integer.
        current_val = int(self.degree_var.get())
        
        # 1. Always update the visual text label instantly
        self.degree_val_label.config(text=str(current_val))

        # 2. THE GATEKEEPER: Only run the heavy math if the integer actually changed
        if current_val != getattr(self, 'last_plotted_degree', -1):
            self.last_plotted_degree = current_val  # Update our tracker
            self.on_plot_click()                    # Fire the backend logic

    def get_axis_limits(self):
        """Safely extracts limits from UI, falling back to defaults if invalid."""
        try: x_min = float(self.x_min_var.get())
        except ValueError: x_min = -1.0
        
        try: x_max = float(self.x_max_var.get())
        except ValueError: x_max = 1.0

        try: y_min = float(self.y_min_var.get())
        except ValueError: y_min = None
        
        try: y_max = float(self.y_max_var.get())
        except ValueError: y_max = None
        
        return x_min, x_max, y_min, y_max

    def get_numeric_coef(self, func_sym, poly_sym, a, b, sigma_sym):
        """Calculates coefficients in milliseconds using SciPy instead of SymPy."""
        # 1. Combine the math symbolically first
        numerator_expr = func_sym * poly_sym * sigma_sym
        denominator_expr = poly_sym * poly_sym * sigma_sym
        
        numerator_expr = powsimp(numerator_expr, combine='exp')
        denominator_expr = powsimp(denominator_expr, combine='exp')

        # 2. Convert to ultra-fast numerical functions
        num_func = lambdify(x, numerator_expr, "numpy")
        den_func = lambdify(x, denominator_expr, "numpy")
        
        # 3. Translate SymPy infinity (oo) to NumPy infinity (np.inf)
        a_num = -np.inf if a == -oo else (np.inf if a == oo else float(a))
        b_num = -np.inf if b == -oo else (np.inf if b == oo else float(b))
        
        # 4. Integrate numerically! (Remember to extract just the value from the tuple)
        num_val, _ = spi.quad(num_func, a_num, b_num)
        den_val, _ = spi.quad(den_func, a_num, b_num)
        
        return num_val / den_val

    def on_plot_click(self):
        """Extracts data and calls the backend logic."""
        mode = self.mode_var.get()
        degree = int(self.degree_var.get())
        
        # --- NEW BASIS SELECTION LOGIC ---
        basis_input_type = self.basis_input_type_var.get()
        
        if basis_input_type == "preset":
            basis_name = self.basis_var.get()
            a, b, sigma = poly_base_params[basis_name]
            cache_key = basis_name # Cache key is just the preset name
        else:
            # Parse the custom inputs! Sympify handles 'oo' natively!
            a = sympify(self.custom_a_var.get())
            b = sympify(self.custom_b_var.get())
            sigma_str = self.custom_sigma_var.get()
            sigma = sympify(sigma_str, locals={'x': x})
            
            basis_name = "Custom Basis" 
            # Create a unique key so the cache knows if you changed the bounds or weight
            cache_key = f"custom_{a}_{b}_{sigma_str}"
        # ---------------------------------

        # ==========================================
        # THE EFFICIENCY BOOST (DICTIONARY CACHING)
        # ==========================================
        # Check if the cache_key exists in our dictionary, and if it has enough terms
        if cache_key not in self.basis_cache or len(self.basis_cache[cache_key]) < degree + 1:
            print(f"Pre-calculating {cache_key} polynomials up to n={self.max_calc_degree}...")
            
            polys = get_ortho_expressions(self.max_calc_degree + 1, a, b, sigma)
            
            with self.cache_lock:
                self.basis_cache[cache_key] = polys
                
            print("Done! Future slider moves will now be instant.")
            
        # Pull the requested number of polynomials from the dictionary
        basis = self.basis_cache[cache_key][:degree + 1]
        # ==========================================

        if mode == "approx":
            input_type = self.input_type_var.get()
            func_str = self.preset_var.get() if input_type == "preset" else self.custom_func_var.get()
            
            # (Remember our sympify fix here!)
            func = sympify(func_str, locals={'x': x}) 
            
            # ==========================================
            # THE COEFFICIENT CACHE
            # ==========================================
            # If the target function OR the basis cache_key changed, wipe the saved coefficients
            if self.cached_func_str != func_str or self.cached_basis_name != cache_key:
                self.cached_coefs = []
                self.cached_func_str = func_str
                self.cached_basis_name = cache_key  # Save the new key!
                
            # Incrementally calculate ONLY the coefficients we are missing!
            while len(self.cached_coefs) <= degree:
                next_n = len(self.cached_coefs)
                
                # Grab the corresponding polynomial from our new DICTIONARY cache
                if next_n >= len(self.basis_cache[cache_key]):
                    print(f"Error, tried to add more coefficients ({next_n}) than there are basis vectors ({len(self.basis_cache[cache_key])})!")
                    break  # Break out to avoid an infinite loop if this happens
                else:
                    poly = self.basis_cache[cache_key][next_n]
                
                    # Calculate and store the single missing coefficient
                    c = self.get_numeric_coef(func, poly, a, b, sigma)
                    self.cached_coefs.append(c)
                
            # Reconstruct the approximation equation instantly using the cached math
            approx_func = sum(c * p for c, p in zip(self.cached_coefs[:degree + 1], basis))
            # ==========================================

            self.draw_plot(func, approx_func, basis_name, degree, a, b)
        else:
            self.draw_basis(basis, basis_name, degree, a, b)

    def draw_plot(self, func, approx_func, basis_name, n, a, b):
        self.ax.clear()
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")
        self.ax.grid(True, linestyle=":", alpha=0.7)

        # ==========================================
        # ADD THE SHADED SAFE ZONE
        # ==========================================
        # Matplotlib crashes on true infinity, so we use a massive finite number instead!
        a_val = -1e9 if a == -oo else float(a)
        b_val = 1e9 if b == oo else float(b)

        self.ax.axvspan(a_val, b_val, color='gold', alpha=0.15, zorder=0, label="Integration Bounds")
        # ==========================================

        # GET CUSTOM LIMITS
        x_min, x_max, y_min, y_max = self.get_axis_limits()

        # USE CUSTOM X LIMITS FOR LINSPACE
        x_vals = np.linspace(x_min, x_max, 2000)
        
        func_numpy = lambdify(x, func, "numpy")
        approx_numpy = lambdify(x, approx_func, "numpy")
        
        y = np.zeros_like(x_vals) + func_numpy(x_vals)
        y_approx = np.zeros_like(x_vals) + approx_numpy(x_vals)
        
        self.ax.set_title(f"Approximation using {basis_name} (n={n})")
        self.ax.plot(x_vals, y, label="Target f(x)", color="black", linewidth=2)
        self.ax.plot(x_vals, y_approx, label=f"P_{n}(x)", color="red", linestyle="--")
        
        # APPLY LIMITS TO MATPLOTLIB
        self.ax.set_xlim(x_min, x_max)
        
        # 1. Calculate the bounds of the target function only
        t_min, t_max = np.min(y), np.max(y)
        
        # 2. Create a 5% visual margin so the line doesn't scrape the edge of the window
        margin = (t_max - t_min) * 0.05
        
        # (Fallback just in case the target is a perfectly flat line, like y=0)
        if margin == 0: 
            margin = 0.5 
            
        # 3. Use custom limits if typed in, otherwise use our calculated target bounds
        final_y_min = y_min if y_min is not None else (t_min - margin)
        final_y_max = y_max if y_max is not None else (t_max + margin)

        self.ax.set_ylim(ymin=final_y_min, ymax=final_y_max)
            
        self.ax.legend()
        self.canvas.draw()

    def draw_basis(self, basis, basis_name, n, a, b):
        self.ax.clear()
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")
        self.ax.grid(True, linestyle=":", alpha=0.7)

        # ==========================================
        # ADD THE SHADED SAFE ZONE
        # ==========================================
        # Matplotlib crashes on true infinity, so we use a massive finite number instead!
        a_val = -1e9 if a == -oo else float(a)
        b_val = 1e9 if b == oo else float(b)

        # axvspan draws a rectangle across the whole Y-axis between two X values
        self.ax.axvspan(a_val, b_val, color='gold', alpha=0.15, zorder=0, label="Integration Bounds")
        # ==========================================

        # GET CUSTOM LIMITS
        x_min, x_max, y_min, y_max = self.get_axis_limits()

        tab20_colors = plt.get_cmap('tab20').colors
        
        # USE CUSTOM X LIMITS FOR LINSPACE
        x_vals = np.linspace(x_min, x_max, 2000)
        
        self.ax.set_title(f"Basis functions {basis_name} (n={n})")
        for bx, basis_func in enumerate(basis):
            basis_func_numpy = lambdify(x, basis_func, "numpy")
            y = np.zeros_like(x_vals) + basis_func_numpy(x_vals)
            self.ax.plot(x_vals, y, label=f"phi_{bx}", color=tab20_colors[(bx*3)%20], linewidth=2)

        # APPLY LIMITS TO MATPLOTLIB
        self.ax.set_xlim(x_min, x_max)
        if y_min is not None or y_max is not None:
            self.ax.set_ylim(ymin=y_min, ymax=y_max)

        self.ax.legend()
        self.canvas.draw()


if __name__ == "__main__":
    app = OrthogonalPolyGUI()
    app.mainloop()