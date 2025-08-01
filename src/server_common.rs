//! Common structures and types used in the client and server
use crate::{IndexHeader, MatchThreshold};
use anyhow::Result;
use reqwest::Client;
use serde::{Deserialize, Serialize};

/// Request structure for filtering minimizers
#[derive(Serialize, Deserialize)]
pub struct FilterRequest {
    /// Prehashed minimizers for input
    pub input: Vec<u64>,

    /// Mininum number (integer) or proportion (float) of minimizer hits for a match
    pub match_threshold: MatchThreshold,
}

/// Response structure for filter results
/// Returns whether this set of minimizers matches the index
#[derive(Serialize, Deserialize)]
pub struct FilterResponse {
    /// Indicates if the input minimizers match the index
    pub index_match: bool,
}

/// Get the header of the index loaded into a remote server
pub async fn get_sever_index_header(server_address: &str) -> Result<IndexHeader> {
    // Create a client to send the minimizers to the server
    let client = Client::new();

    // Send the minimizers as a POST request
    let response = client
        .get(server_address.to_owned() + "/index_header")
        .send()
        .await?;

    // Check if the response indicates a match
    if response.status().is_success() {
        Ok(response.json::<IndexHeader>().await?)
    } else {
        Err(anyhow::anyhow!(
            "Server returned an error: {}",
            response.status()
        ))
    }
}
