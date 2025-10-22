# Copyright 2025 Matias Senger. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from pathlib import Path
import subprocess
from logging import getLogger, basicConfig, INFO
import pandas
import plotly.express as px

logger = getLogger(__file__)

def run_palace(file:Path,n_processors:int=1):
    cmd = ['palace','-np',str(n_processors),str(file)]
    logger.info(f'Running command `{" ".join(cmd)}`...')
    result = subprocess.run(
        cmd, # This is just the command, e.g. `palace -np 1 /path/to/coaxial_matched.json`.
        check = True, # Check that there were no errors.
        stdout = subprocess.DEVNULL, # Hide prints from `palace...`.
        stderr = subprocess.DEVNULL, # Hide prints from `palace...`.
    )

def run_simulations():
    # Run each of the three simulations:
    for fname in {'coaxial_matched.json','coaxial_open.json','coaxial_short.json'}:
        run_palace(Path(__file__).parent/fname)

def plot_results():
    logger.info('Plotting results...')
    PLOTS_DIR = Path(__file__).parent/'postpro/plots'
    PLOTS_DIR.mkdir(exist_ok=True)
    # Read the data:
    data = []
    for ending in {'matched','open','short'}:
        __ = []
        for fname in {'domain-E.csv','port-I.csv','port-V.csv'}:
            _ = pandas.read_csv(Path(__file__).parent/'postpro'/ending/fname, skipinitialspace=True)
            _['coax ending'] = ending
            _.set_index(['t (ns)','coax ending'], inplace=True)
            __.append(_)
        __ = pandas.concat(__, axis=1)
        data.append(__)
    data = pandas.concat(data)
    # Do the plots:
    for variable in {'V[1] (V)','I[1] (A)','E_elec (J)','E_mag (J)'}:
        fig = px.line(
            title = f'{variable.split(" ")[0]} vs time',
            data_frame = data.sort_index().reset_index(drop=False),
            x = 't (ns)',
            y = variable,
            facet_row = 'coax ending',
        )
        fig.write_html(PLOTS_DIR/f'{variable}.html')
    logger.info(f'Plots available in {PLOTS_DIR}')

def main():
    run_simulations() # Here we call *palace*.
    plot_results() # Do some plots to visualize the results.

if __name__ == '__main__':
    # Configure the logger:
    basicConfig(
        level = INFO,
        format = '%(asctime)s %(levelname)s - %(message)s',
        datefmt = '%b-%d %H:%M:%S',
    )
    # Run the interesting thing:
    main()
