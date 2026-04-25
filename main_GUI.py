import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

class OrthogonalPolyGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Orthogonal Polynomial Visualizer")
        self.geometry("900x700")
        self.minsize(700, 500)

        # Configure grid to allow resizing
        self.columnconfigure(0, weight=1) # Control panel
        self.columnconfigure(1, weight=4) # Plot area
        self.rowconfigure(0, weight=1)

        self.create_widgets()
        self.update_ui_state() # Initialize the correct disabled/enabled states

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

        # 2. Dropdown Menu for Polynomial Basis
        ttk.Label(control_frame, text="Polynomial Basis:").pack(anchor="w", pady=(0, 5))
        self.basis_var = tk.StringVar(value="Legendre")
        self.basis_dropdown = ttk.Combobox(control_frame, textvariable=self.basis_var, state="readonly")
        self.basis_dropdown['values'] = ("Legendre", "Hermite", "Chebyshev", "Laguerre")
        self.basis_dropdown.pack(fill="x", pady=(0, 15))

        # 3. Function Input Section
        self.func_frame = ttk.LabelFrame(control_frame, text="Target Function", padding="10")
        self.func_frame.pack(fill="x", pady=(0, 15))

        self.input_type_var = tk.StringVar(value="preset")
        self.input_type_var.trace_add("write", self.update_ui_state)

        # 3a. Preset Functions
        ttk.Radiobutton(self.func_frame, text="Preset Function", variable=self.input_type_var, value="preset").pack(anchor="w")
        self.preset_var = tk.StringVar(value="sin(x)")
        self.preset_dropdown = ttk.Combobox(self.func_frame, textvariable=self.preset_var, state="readonly")
        self.preset_dropdown['values'] = ("sin(x)", "cos(x)", "e^x", "Step Function", "Absolute Value |x|")
        self.preset_dropdown.pack(fill="x", padx=20, pady=(0, 10))

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
            slider_frame, from_=1, to=15, orient="horizontal",
            variable=self.degree_var, command=self.update_slider_label
        )
        self.degree_slider.pack(side="left", fill="x", expand=True)
        
        self.degree_val_label = ttk.Label(slider_frame, text="5")
        self.degree_val_label.pack(side="right", padx=(10, 0))

        # 5. Update Button
        self.update_btn = ttk.Button(control_frame, text="Update Plot", command=self.on_plot_click)
        self.update_btn.pack(fill="x", pady=(10, 0))

        # --- Plot Area (Right Side) ---
        plot_frame = ttk.Frame(self, padding="10")
        plot_frame.grid(row=0, column=1, sticky="nsew")

        self.figure = Figure(figsize=(6, 5), dpi=100)
        self.ax = self.figure.add_subplot(111)
        
        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        self.draw_plot()

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

    def update_slider_label(self, event):
        val = int(self.degree_var.get())
        self.degree_val_label.config(text=str(val))

    def on_plot_click(self):
        """Extracts data and calls the backend logic."""
        mode = self.mode_var.get()
        basis = self.basis_var.get()
        degree = int(self.degree_var.get())
        
        if mode == "approx":
            input_type = self.input_type_var.get()
            func_str = self.preset_var.get() if input_type == "preset" else self.custom_func_var.get()
            print(f"BACKEND CALL -> Mode: {mode}, Basis: {basis}, Func: {func_str}, Degree: {degree}")
        else:
            func_str = None
            print(f"BACKEND CALL -> Mode: {mode}, Basis: {basis}, Displaying first {degree} polynomials")

        # ==========================================
        # BACKEND INTEGRATION POINT
        # ==========================================
        # if mode == "approx":
        #     x, y_true, y_approx = backend.approximate(basis, func_str, degree)
        #     self.draw_plot(mode, x, y_true=y_true, y_approx=y_approx)
        # elif mode == "basis":
        #     x, list_of_y_arrays = backend.get_basis_polynomials(basis, degree)
        #     self.draw_plot(mode, x, basis_ys=list_of_y_arrays)
        # ==========================================

        self.draw_plot()

    def draw_plot(self, mode=None, x=None, y_true=None, y_approx=None, basis_ys=None):
        """Clears the old plot and draws the new fake data based on mode."""
        self.ax.clear()
        
        if mode is None:
            mode = self.mode_var.get()
        
        basis = self.basis_var.get()
        degree = int(self.degree_var.get())
        
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")
        self.ax.grid(True, linestyle=":", alpha=0.7)

        if x is None:
            x = np.linspace(-1, 1, 200)

        if mode == "approx":
            self.ax.set_title(f"Approximation using {basis} (n={degree})")
            
            # Fake Data for Approximation
            if y_true is None: y_true = np.sin(x * np.pi)
            if y_approx is None: y_approx = y_true + np.cos(x * degree * np.pi) * (1.0 / degree)

            self.ax.plot(x, y_true, label="Target f(x)", color="black", linewidth=2)
            self.ax.plot(x, y_approx, label=f"P_{degree}(x)", color="red", linestyle="--")
            self.ax.legend()

        elif mode == "basis":
            self.ax.set_title(f"First {degree} {basis} Polynomials")
            
            # Fake Data for Basis Polynomials (Just plotting x^0, x^1, x^2...)
            if basis_ys is None:
                basis_ys = [x**i for i in range(degree)]
                
            for i, y_vals in enumerate(basis_ys):
                self.ax.plot(x, y_vals, label=f"n={i}")
            
            # Only show legend if there aren't too many lines
            if degree <= 10:
                self.ax.legend(loc="upper right", bbox_to_anchor=(1.15, 1))

        self.canvas.draw()

if __name__ == "__main__":
    app = OrthogonalPolyGUI()
    app.mainloop()