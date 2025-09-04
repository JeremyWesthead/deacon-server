//! Common structures and types used in the client and server
use crate::IndexHeader;
use anyhow::Result;
use reqwest::blocking::Client;
use serde::{Deserialize, Serialize};

/// Request structure for filtering minimizers
#[derive(Serialize, Deserialize, Clone)]
pub struct FilterRequest {
    /// Prehashed minimizers for input
    pub hashes: Vec<u64>,
    /// Positions of the minimizers in the input sequence
    pub valid_positions: Vec<u32>,
    /// Effective sequence of the input (e.g., nucleotide sequence)
    pub effective_seqs: Vec<Vec<u8>>,
    /// Length of the k-mer used for minimization
    pub kmer_length: usize,
    /// Whether running in debug mode
    pub debug: bool,
}

/// Response structure for filter results
/// Returns whether this set of minimizers should be output
#[derive(Serialize, Deserialize, Clone)]
pub struct FilterResponse {
    /// The k-mers that matched the index
    pub hit_kmers: Vec<String>,
    /// Number of index hits
    pub hit_count: usize,
}

/// Get the header of the index loaded into a remote server
/// Required in order to ensure that the locally computed minimizers match
/// the kmer length and window size
pub fn get_server_index_header(server_address: &str) -> Result<IndexHeader> {
    // Create a client to send the minimizers to the server
    let client = Client::new();

    // Send the minimizers as a POST request
    let response = client
        .get(server_address.to_owned() + "/index_header")
        .send()?;

    // Check if the response indicates a match
    if response.status().is_success() {
        Ok(response.json::<IndexHeader>()?)
    } else {
        Err(anyhow::anyhow!(
            "Server returned an error: {}",
            response.status()
        ))
    }
}
