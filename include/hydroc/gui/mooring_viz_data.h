/// @file mooring_viz_data.h
/// @brief Lightweight data types for mooring-line visualization.
///
/// Header-only with no MoorDyn or Eigen dependencies so that both the
/// physics library (HydroChrono) and the visualization library
/// (HydroChronoGUI) can include it without pulling in heavy headers.
#pragma once

#include <array>
#include <functional>
#include <vector>

namespace hydroc {
namespace gui {

/// MoorDyn connection-point classification.
/// Values intentionally match the MoorDyn C API convention.
enum class MooringPointType : int {
    kCoupled = -1,  ///< Vessel fairlead (driven by the coupled body).
    kFree    =  0,  ///< Free point / clump weight / buoy.
    kFixed   =  1,  ///< Seabed anchor or fixed attachment.
};

/// Sentinel value used for intermediate (non-endpoint) node markers.
inline constexpr int kIntermediateNode = -99;

/// Per-line visualization payload: ordered node positions along the line
/// plus the MoorDyn point type at each endpoint.
struct MooringLineVizData {
    std::vector<std::array<double, 3>> node_positions;
    std::vector<double> node_tensions;      ///< Per-node tension magnitude [N].
    double line_max_tension = 0.0;          ///< Max tension magnitude along line [N].
    MooringPointType start_point_type = MooringPointType::kFree;
    MooringPointType end_point_type   = MooringPointType::kFree;
};

/// Callback the GUI invokes each frame to fetch current line geometry.
using MooringVizProvider =
    std::function<std::vector<MooringLineVizData>()>;

}  // namespace gui
}  // namespace hydroc
