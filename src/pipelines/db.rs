use std::path::PathBuf;
use std::sync::Arc;
use std::fs;
use anyhow::anyhow;
use log::info;
use crate::config::defs::{PipelineError, RunConfig};
use crate::utils::file::{file_path_manipulator, validate_file_inputs};
use crate::utils::taxonomy::{build_taxid_lineages_db, build_accession2taxid_db, build_fst_acc2taxid};

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

    let db_out_path = file_path_manipulator(&PathBuf::from("taxid_lineage"), Some(&config.out_dir), None, Some(".bincode"), "_");

    build_taxid_lineages_db(&nodes_dmp_path, &merged_dmp_path, &db_out_path)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    Ok(())
}


pub async fn accession2taxid_db(config: Arc<RunConfig>) -> anyhow::Result<(), PipelineError> {
    let cwd = std::env::current_dir().map_err(|e| PipelineError::Other(e.into()))?;

    info!("Building accession2taxid db.");

    let (mandatory_path, _file2_path, _no_ext_sample_base_buf, _no_ext_sample_base) =
        validate_file_inputs(&config, &cwd)?;

    let gz_paths: Vec<PathBuf> = if mandatory_path.is_dir() {
        // if dir, all files inside assumed to be acc2tax source files
        let mut entries = tokio::fs::read_dir(&mandatory_path).await
            .map_err(|e| PipelineError::Other(anyhow!("cannot read directory {}: {}", mandatory_path.display(), e)))?;

        let mut paths = Vec::new();
        while let Some(entry) = entries.next_entry().await? {
            if entry.file_type().await?.is_file() {
                paths.push(entry.path());
            }
        }
        info!("Found {} source files in directory {}", paths.len(), mandatory_path.display());
        paths
    } else if mandatory_path.is_file() {
        info!("Using single source file {}", mandatory_path.display());
        vec![mandatory_path.clone()]
    } else {
        return Err(PipelineError::Other(anyhow!(
            "mandatory input must be a file or directory: {}",
            mandatory_path.display()
        )));
    };

    if gz_paths.is_empty() {
        return Err(PipelineError::Other(anyhow!(
            "no accession2taxid source files found in {}",
            mandatory_path.display()
        )));
    }

    let nt_path: Option<PathBuf> = config.args.nt.clone().map(PathBuf::from);
    let nr_path: Option<PathBuf> = config.args.nr.clone().map(PathBuf::from);

    let db_out_path = file_path_manipulator(&PathBuf::from("accession2_taxid"), Some(&cwd), None, Some(".fst"), "_");
    eprintln!("Writing to DB {:?}", db_out_path);

    // build_accession2taxid_db(&gz_paths, nt_path.as_ref(), nr_path.as_ref(), &db_out_path)
    //     .await
    //     .map_err(|e| PipelineError::Other(e.into()))?;

    build_fst_acc2taxid(&gz_paths, nt_path.as_ref(), nr_path.as_ref(), &db_out_path)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;



    Ok(())
}