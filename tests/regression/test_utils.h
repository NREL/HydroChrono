#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include <chrono/utils/ChValidation.h>
#include <chrono_postprocess/ChGnuPlot.h>

void PlotValidation(const std::string& out_script,
                    const std::string& title,
                    const chrono::utils::ChValidation::Headers& headers,
                    chrono::utils::ChValidation::Data& ref_data,
                    chrono::utils::ChValidation::Data& res_data,
                    double simulationDuration) {
    int n_plots = (int)res_data.size() - 1;

    chrono::postprocess::ChGnuPlot gplot(out_script);
    gplot.SetCanvasSize(1200, 500 * n_plots);
    gplot.SetTitle("{/:Bold " + title + "}");
    gplot.SetCommand("set multiplot layout " + std::to_string(n_plots) + ", 1");

    for (int i = 1; i <= n_plots; i++) {
        gplot.SetRangeX(0, simulationDuration);
        gplot.SetCommand("set xlabel '" + headers[0] + "' noenhanced");

        gplot.SetCommand("set ylabel '" + headers[i] + "' noenhanced");        
        gplot.Plot(ref_data[0], ref_data[i], "Reference", "with lines lt rgb '#0055FF' lw 4");
        gplot.Plot(res_data[0], res_data[i], "Simulation", "with lines lt rgb '#FF3300' lw 4 dashtype 2");

        auto n_err                      = std::min(res_data[0].size(), ref_data[0].size());
        chrono::ChVectorDynamic<> t_err = ref_data[0].head(n_err);
        chrono::ChVectorDynamic<> err   = ref_data[i].head(n_err) - res_data[i].head(n_err);
        if (err.norm() < 1e-4) {
            gplot.SetRangeY2(-0.01, +0.01);
            gplot.SetCommand("set y2tics -0.01 , 0.01, 0.01");
        }
        gplot.SetLabelY2("error");
        gplot.Plot(t_err, err, "Error (ref - sim)", "axes x1y2 with lines lt rgb '#44CC44' lw 4");

        gplot.SetLegend("box lt -1 lw 1");
        gplot.SetGrid();

        gplot.FlushPlots();
    }
}

#endif
