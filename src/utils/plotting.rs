use anyhow::Result;
use plotters::prelude::*;
use std::collections::HashMap;
use std::path::Path;

pub fn plot_depths(depth_map: &HashMap<usize, u32>, sample_name: &str, output_path: &Path) -> Result<()> {
    // Convert depth_map to Vec, ensuring all positions are included
    let max_pos = depth_map.keys().max().copied().unwrap_or(1);
    let depths: Vec<u32> = (1..=max_pos).map(|i| *depth_map.get(&i).unwrap_or(&0)).collect();

    if depths.is_empty() {
        return Err(anyhow::anyhow!("No depth data available for plotting"));
    }

    let root = BitMapBackend::new(output_path, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    // Find max
    let max_depth = depths.iter().max().copied().unwrap_or(1).max(1);

    let mut chart = ChartBuilder::on(&root)
        .caption(sample_name, ("sans-serif", 20))
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(1..max_pos, LogRange(1u32..max_depth + 1))?;

    chart.configure_mesh()
        .x_desc("position")
        .y_desc("depth")
        .draw()?;

    chart.draw_series(LineSeries::new(
        depths.iter().enumerate().map(|(i, &d)| (i + 1, d)),
        &BLUE,
    ))?;

    root.present()?;
    Ok(())
}