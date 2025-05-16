# IntestinalIntussusceptionApp\_e

A MATLAB-based computational model for simulating pediatric small intestine peristalsis and elasticity to predict intussusception risk under post-gastroenteritis conditions.

## Features

* **1D Fluid–Structure Interaction**: Discretizes a 300 cm small intestine into 100 segments, modeling viscoelastic wall mechanics and luminal content propagation.
* **Customizable Parameters**: Adjust simulation time, peristaltic speed profiles, elasticity spectra, and risk thresholds via a user-friendly GUI.
* **Seven Speed Patterns**: Linear, segmented, physiological, linear/sigmoid/exponential transitions, and step decline near the terminal ileum.
* **Elasticity Distributions**: Uniform, linear gradient, sigmoid transition, anatomical regions, and terminal-ileum focus.
* **Risk Index Calculation**: Combines diameter ratio, speed gradient, and elasticity factors to flag potential intussusception events.

## Requirements

* MATLAB R2022a or later
* Signal Processing Toolbox
* Optimization Toolbox (for parameter sweeps)

## Installation

1. Clone this repository:

   ```sh
   git clone https://github.com/yourusername/IntestinalIntussusceptionModel.git
   ```
2. Open MATLAB and navigate to the project folder:

   ```matlab
   cd('path/to/IntestinalIntussusceptionModel')
   IntestinalIntussusceptionApp_e.mlapp
   ```
3. Ensure all dependencies are installed via MATLAB Add-On Explorer.

## Usage

1. Launch the App by double-clicking `IntestinalIntussusceptionApp_e.mlapp` or running:

   ```matlab
   app = IntestinalIntussusceptionApp_e;
   ```
2. Configure parameters in the **Settings** panel:

   * Simulation Time, Onset Time
   * Speed Profiles & Variability
   * Elasticity Spectra & Variability
   * Risk Thresholds
3. Click **Run Simulation**.
4. View results in the **Visualization** tab:

   * Speed & Elasticity profiles
   * Risk heatmaps along the intestinal length
   * Event log for detected intussusception segments
5. Export results as CSV or PNG via the **Export** menu.

## Examples

Pre-configured examples are provided in the `examples/` directory:

* `example_speed_profiles.m` - Demonstrates all seven speed patterns.
* `example_sensitivity_analysis.m` - Shows elasticity gradient vs. value factor sensitivity.

## Contributing

Contributions are welcome! Please fork the repository and submit a pull request. For major changes, open an issue first to discuss proposed improvements.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Citation

If you use this code in your research, please cite:

> Hsieh, M.-Y., et al. Predicting Pediatric Intussusception Risk via Computational Modeling of Intestinal Elasticity and Peristaltic Speed. *Journal Name*, Volume(Issue)\:Pages (2025).
