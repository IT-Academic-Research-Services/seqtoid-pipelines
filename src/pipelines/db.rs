use std::path::PathBuf;
use std::sync::Arc;
use std::fs;
use anyhow::anyhow;
use log::info;
use crate::config::defs::{PipelineError, RunConfig};
use crate::utils::file::{file_path_manipulator, validate_file_inputs};
use crate::utils::taxonomy::{build_taxid_lineages_db, build_accession2taxid_db, unpack_lineage_bytes};

pub async fn taxid_lineages_db(config: Arc<RunConfig>) -> anyhow::Result<(), PipelineError> {
    let cwd = std::env::current_dir().map_err(|e| PipelineError::Other(e.into()))?;

    info!("Building taxid lineages db.");

    let (taxid_dir_path, _file2_path, no_ext_sample_base_buf, _no_ext_sample_base) = validate_file_inputs(&config, &cwd)?;
    eprintln!("base {}", no_ext_sample_base_buf.display());
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

    let db_out_path = file_path_manipulator(&PathBuf::from("taxid_lineage"), Some(&config.out_dir), None, Some("sled.db"), "_");

    build_taxid_lineages_db(&nodes_dmp_path, &merged_dmp_path, &db_out_path)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    if config.args.verbose {
        let db = sled::open(db_out_path).map_err(|e| PipelineError::Other(e.into()))?;
        let tree = db.open_tree("lineages").map_err(|e| PipelineError::Other(e.into()))?;
        let count = tree.iter().count();
        println!("Total entries: {}", count);

        for taxid in config.args.test_taxids.clone() {
            if let Some(ivec) = tree.get(&taxid.to_le_bytes()).map_err(|e| PipelineError::Other(e.into()))? {
                let lineage = unpack_lineage_bytes(&ivec)?;
                println!("Taxid {} lineage: {:?}", taxid, lineage);
            } else {
                println!("Taxid {} not found", taxid);
            }
        }

    }

    Ok(())
}


pub async fn accession2taxid_db(config: Arc<RunConfig>) -> anyhow::Result<(), PipelineError> {
    let cwd = std::env::current_dir().map_err(|e| PipelineError::Other(e.into()))?;

    info!("Building accession2taxid db.");

    let (accession2taxid_file_path, _file2_path, no_ext_sample_base_buf, _no_ext_sample_base) = validate_file_inputs(&config, &cwd)?;

    let db_out_path = file_path_manipulator(&PathBuf::from("accession2_taxid"), Some(&cwd), None, Some("sled.db"), "_");
    eprintln!("Writing to DB {:?}", db_out_path);

    build_accession2taxid_db(&accession2taxid_file_path, &db_out_path)
        .await
        .map_err(|e| PipelineError::Other(e.into()))?;

    if config.args.verbose {
        let db = sled::open(db_out_path).map_err(|e| PipelineError::Other(e.into()))?;
        let tree = db.open_tree("acc2taxid")?;
        let count = tree.iter().count();
        println!("Total entries: {}", count);

        for acc in config.args.test_accessions {
            if let Some(ivec) = tree.get(acc.as_bytes())? {
                let taxid = i32::from_le_bytes(ivec[..4].try_into().map_err(|_| anyhow!("corrupt taxid"))?);
                println!("Accession {} taxid: {}", acc, taxid);
            } else {
                println!("Accession {} not found", acc);
            }
        }
    }

    Ok(())
}