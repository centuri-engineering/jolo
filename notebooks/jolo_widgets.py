import os
from pathlib import Path
from threading import Thread, Event
from multiprocessing import Process
import signal

import numpy as np
import pandas as pd
import altair as alt
import ipywidgets as widgets

from pygtftk.plugins.ologram import ologram
from pygtftk.plugins.ologram import make_parser
import pygtftk.utils


def feature_graph(tsv_file):
    """This could use some polishing
    """
    result = pd.read_csv(tsv_file, sep="\t")
    fc = result["nb_intersections_true"] / (
        result["nb_intersections_expectation_shuffled"] + 1
    )

    result["signif_color"] = "#b3b3b3"

    result.loc[
        (result.nb_intersections_pvalue < 0.05) & (fc < 1), "signif_color"
    ] = "#ffa64d"
    result.loc[
        (result.nb_intersections_pvalue < 0.05) & (fc > 1), "signif_color"
    ] = "#6cc67b"
    result.loc[
        (result.nb_intersections_pvalue < 1e-10) & (fc < 1), "signif_color"
    ] = "#cc6600"
    result.loc[
        (result.nb_intersections_pvalue < 1e-10) & (fc > 1), "signif_color"
    ] = "#3c9040"

    result["err_min"] = (
        result["nb_intersections_expectation_shuffled"]
        - result["nb_intersections_variance_shuffled"] ** 0.5
    )

    result["err_max"] = (
        result["nb_intersections_expectation_shuffled"]
        + result["nb_intersections_variance_shuffled"] ** 0.5
    )

    chart = alt.Chart(
        result, width=600, height=400, title="Total overlap length per region type"
    )

    x = alt.X("feature_type:O", title="Feature type")
    y = alt.Y(
        "nb_intersections_expectation_shuffled", title="Nb of overlapping base pairs"
    )

    bar_width = 15
    bars = chart.mark_bar(color="#cccccc", width=bar_width,).encode(
        x=x, y=y, tooltip="nb_intersections_pvalue"
    )

    errs = chart.mark_rule().encode(x=x, y="err_min", y2="err_max",)

    y_true = alt.Y("nb_intersections_true")
    barrs2 = chart.mark_bar(color="#67b7e3", width=bar_width, xOffset=bar_width).encode(
        x=x, y=y_true, tooltip="nb_intersections_pvalue"
    )

    pvals = chart.mark_text(yOffset=-20).encode(
        x=x, y=y_true, text="nb_intersections_pvalue",
    )

    return bars + errs + barrs2 + pvals


# TODO: add a selector for feature sorting:

# sort_features_opts = [
#     None,
#     "nb_intersections_expectation_shuffled",
#     "nb_intersections_variance_shuffled",
#     "nb_intersections_negbinom_fit_quality",
#     "nb_intersections_log2_fold_change",
#     "nb_intersections_true",
#     "nb_intersections_pvalue",
#     "summed_bp_overlaps_expectation_shuffled",
#     "summed_bp_overlaps_variance_shuffled",
#     "summed_bp_overlaps_negbinom_fit_quality",
#     "summed_bp_overlaps_log2_fold_change",
#     "summed_bp_overlaps_true",
#     "summed_bp_overlaps_pvalue",
# ]


def _browse_args(parser):
    """relies on argparse.ArgumentParser internals
    to gather arguments and help strings
    (note that one should not in theory access private members
    but here it simplifies a lot of things)
    """
    all_args = {}
    for action in parser._actions:
        key = action.option_strings[1].lstrip("--")
        short = action.option_strings[0]
        all_args[key] = (short, action.help)
    return all_args


class BedList(widgets.Accordion):
    def __init__(self, beds, **kwargs):
        self.beds = beds
        self.check_boxes = [
            widgets.Checkbox(value=False, description=bed[0], disabled=False)
            for bed in beds
        ]

        children = [widgets.VBox(self.check_boxes, **kwargs)]
        super().__init__(children=children, selected_index=None)
        super().set_title(
            0, "Select more beds",
        )

    @property
    def selected(self):
        return [bed[1] for bed, box in zip(self.beds, self.check_boxes) if box.value]

    @property
    def value(self):
        return self.selected


class OlogramCmd(widgets.VBox):
    def __init__(self, input_dir, output_dir):

        style = {"description_width": "200px"}
        self.output_dir = output_dir
        self.parser = make_parser()
        self.all_args = _browse_args(self.parser)

        self.gtfs = [(f.name, f.as_posix()) for f in Path(input_dir).glob("*.gtf*")]
        self.beds = [(f.name, f.as_posix()) for f in Path(input_dir).glob("*.bed")]
        self.genomes = [
            (f.name, f.as_posix()) for f in Path(input_dir).glob("*.genome")
        ]

        self.arg_widgets = {
            "inputfile": widgets.Dropdown(
                options=self.gtfs, description="Input File", style=style,
            ),
            "peak-file": widgets.Dropdown(
                options=self.beds, description="Peak File", style=style
            ),
            "chrom-info": widgets.Dropdown(
                options=self.genomes, description="Chromosome info", style=style
            ),
            "upstream": widgets.BoundedIntText(
                value=1500,
                step=100,
                min=0,
                max=10000,  # set reasonable values here
                description="Upstream extension",
                style=style,
            ),
            "downstream": widgets.BoundedIntText(
                value=1500,
                step=100,
                min=0,
                max=10000,  # TODO set reasonable values here
                description="Downstream extension",
                style=style,
            ),
            # Maybe this should not be configurable
            "nb-threads": widgets.BoundedIntText(
                value=1,
                min=1,
                max=os.cpu_count(),
                description="Number of threads to use",
                style=style,
            ),
            "more-bed": BedList(
                self.beds, description="Select more beds", style=style,
            ),
        }
        self.cmd_box = widgets.HTML(
            value=f"""<h4>Generated command:\n\n</h4><p><code>{self.cmd}<code/></h4>"""
        )

        for k, w in self.arg_widgets.items():
            w.description_tooltip = self.all_args.get(k, ["", ""])[1]
            w.observe(self.on_cmd_change)

        options = widgets.VBox(list(self.arg_widgets.values()))
        self.run_btn = widgets.Button(description="run")
        self.run_btn.button_style = "info"
        self.run_btn.on_click(self.run)
        self.interupt_button = widgets.Button(description="stop")
        self.interupt_button.on_click(self.stop)
        self.output = widgets.Output()
        self.current_process = Process(target=lambda: print("Not set"))
        super().__init__(
            [
                options,
                widgets.HBox([self.run_btn, self.interupt_button]),
                self.cmd_box,
                self.output,
            ],
        )

    def on_cmd_change(self, change):

        value = f"""<h4>Generated command:\n\n</h4><p><code>{self.cmd}</code></p>"""
        self.cmd_box.value = value
        self.run_btn.description = "run"
        self.run_btn.button_style = "info"

    @property
    def args(self):
        _args = []
        for k, w in self.arg_widgets.items():
            _args.append(f"--{k}")
            if isinstance(w.value, list):
                _args.extend(w.value)
            else:
                _args.append(str(w.value))
        _args.append("-x")
        _args.append(f"--outputdir")
        _args.append(f"{self.output_dir}")
        return _args

    @property
    def kwargs(self):
        try:
            return self.parser.parse_args(self.args).__dict__
        except SystemExit:
            print("Invalid Command, parser tried to exit")
            return None
        except Exception as e:
            print("Invalid command")
            print(e)
            return None

    @property
    def cmd(self):
        _cmd = ["gtftk", "ologram"]
        for k, w in self.arg_widgets.items():
            if isinstance(w.value, bool) and w.value:
                _cmd.append(self.all_args.get(k, [""])[0])
            else:
                # check for empty lists
                if not w.value:
                    continue
                _cmd.append(self.all_args.get(k, [""])[0])
                if isinstance(w.value, list):
                    _cmd.append(" ".join(w.value))
                else:
                    _cmd.append(str(w.value))
        return " ".join(_cmd)

    def run(self, btn):

        if self.current_process.is_alive():
            return

        self.current_process.close()
        with self.output:
            if not self.kwargs:
                self.run_btn.button_style = "danger"
                self.run_btn.description = "errored"
                return
            self.run_btn.description = "running ..."
            self.run_btn.button_style = "warning"
            self.current_process = Process(target=self._call_ologram, name="main")
            self.current_process.start()
            self.output.clear_output()
            return 0

    def _call_ologram(self):
        ologram(**self.kwargs)
        self.run_btn.description = "done"
        self.run_btn.button_style = "success"

    def stop(self, btn):
        if not self.current_process.is_alive():
            return

        self.run_btn.description = "stopping"
        self.current_process.terminate()
        self.current_process.join()
        self.current_process.close()
        self.run_btn.description = "stopped"
        self.run_btn.button_style = "danger"
