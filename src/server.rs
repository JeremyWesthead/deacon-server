//! Functionality to create a server endpoint which can be used to filter based on a pre-loaded index
use std::path::PathBuf;
use std::sync::{OnceLock, RwLock};

use crate::filter::get_index_matches;
use crate::index::{IndexHeader, load_minimizer_hashes};
use crate::server_common::{FilterRequest, FilterResponse};
use axum::{
    Json, Router,
    extract::DefaultBodyLimit,
    routing::{get, post},
};
use rustc_hash::FxHashSet;

/// Shared index file between endpoint calls.
static INDEX: RwLock<Option<FxHashSet<u64>>> = RwLock::new(None);

/// Shared index header between endpoint calls.
static INDEX_HEADER: RwLock<Option<IndexHeader>> = RwLock::new(None);

/// Just for ensuring we get a single tracing setup.
/// Mostly needed as tests otherwise try to spawn multiple
static TRACING: OnceLock<()> = OnceLock::new();

/// Starts the server with the given index path and port.
/// To log the server's connections, set `RUST_LOG=trace` in your environment variables.
pub async fn run_server(index_path: PathBuf, port: u16) {
    // initialize tracing
    TRACING.get_or_init(|| {
        tracing_subscriber::fmt::init();
    });

    eprintln!("Loading index from: {}", index_path.display());
    // Load the index before starting the server to ensure it's available for requests
    load_index(index_path);

    let app = Router::new()
        .route("/", get(root))
        .route("/index_header", get(index_header))
        .route("/get_index_matches", post(index_matches))
        // Increase the body limit to 2GB to ensure we don't error on large payloads
        .layer(DefaultBodyLimit::max(2147483648));

    // run our app with hyper, listening globally
    let listener = tokio::net::TcpListener::bind("0.0.0.0:".to_owned() + &port.to_string())
        .await
        .unwrap();
    axum::serve(listener, app).await.unwrap();
}

fn load_index(index_path: PathBuf) {
    let (index, header) = load_minimizer_hashes(&Some(&index_path), &None).expect("Failed to load index");
    match index {
        Some(index) => {
            INDEX.write().unwrap().replace(index);
            INDEX_HEADER.write().unwrap().replace(header);
        },
        None => panic!("Failed to load index from {index_path:?}"),
    }
}


pub async fn root() -> String {
    let index_header = INDEX_HEADER.read().unwrap();
    let index = INDEX.read().unwrap();
    match (&*index_header, &*index) {
        (Some(header), Some(index)) => format!("Index loaded with {} minimizers, and header: {:?}", index.len(), header),
        _ => "No index loaded".to_string(),
    }
} 

pub async fn index_header() -> Json<IndexHeader> {
    let index_header = INDEX_HEADER.read().unwrap();
    match &*index_header {
        Some(header) => Json(header.clone()),
        None => panic!("No index header loaded"),
    }
}

pub async fn index_matches(Json(request): Json<FilterRequest>,) -> Json<FilterResponse> {
    let index = INDEX.read().unwrap().clone().unwrap();
    let (hit_kmers, hit_count) = get_index_matches(&index, &request.hashes, request.valid_positions, request.effective_seqs, request.kmer_length, request.debug);
    Json(FilterResponse {
        hit_kmers,
        hit_count,
    })
} 


