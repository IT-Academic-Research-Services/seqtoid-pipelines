use std::path::PathBuf;
use std::sync::Arc;
use std::fs;
use anyhow::anyhow;
use log::info;
use crate::config::defs::{PipelineError, RunConfig};
use crate::utils::file::{file_path_manipulator, validate_file_inputs};
use crate::utils::taxonomy::{build_taxid_lineages_db, build_accession2taxid_db};

pub async fn taxid_lineages_db(config: Arc<RunConfig>) -> anyhow::Result<(), PipelineError> {
    let cwd = std::env::current_dir().map_err(|e| PipelineError::Other(e.into()))?;

    info!("Building taxid lineages db.");

    let (taxid_dir_path, _file2_path, no_ext_sample_base_buf, _no_ext_sample_base) = validate_file_inputs(&config, &cwd)?;
    let metadata = fs::metadata(taxid_dir_path.clone())?;
    let file_type = metadata.file_type();
    if !file_type.is_dir() {
        return Err(PipelineError::NotDirectory(taxid_dir_path));
    }

    let nodes_dmp_path = taxid_dir_path.join("nodes.dmp");
    if !nodes_dmp_path.exists() {
        return Err(PipelineError::FileNotFound(nodes_dmp_path));
    }

    let merged_dmp_path = taxid_dir_path.join("merged.dmp");
    if !merged_dmp_path.exists() {
        return Err(PipelineError::FileNotFound(nodes_dmp_path));
    }

    let db_out_path = file_path_manipulator(&no_ext_sample_base_buf, Some(&cwd), None, Some("taxid_lineage_sled.db"), "_");

    eprintln!("Writing to DB {:?}", db_out_path);
    build_taxid_lineages_db(&nodes_dmp_path, &merged_dmp_path, &db_out_path)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    Ok(())
}


pub async fn accession2taxid_db(config: Arc<RunConfig>) -> anyhow::Result<(), PipelineError> {
    let cwd = std::env::current_dir().map_err(|e| PipelineError::Other(e.into()))?;

    info!("Building accession2taxid db.");

    let (accession2taxid_file_path, _file2_path, no_ext_sample_base_buf, _no_ext_sample_base) = validate_file_inputs(&config, &cwd)?;

    let db_out_path = file_path_manipulator(&no_ext_sample_base_buf, Some(&cwd), None, Some("accession2_taxid_sled.db"), "_");
    eprintln!("Writing to DB {:?}", db_out_path);

    build_accession2taxid_db(&accession2taxid_file_path, &db_out_path)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    Ok(())
}