//! Functionality to create a server endpoint which can be used to filter based on a pre-loaded index
use std::path::PathBuf;
use std::sync::{Mutex, OnceLock};

use crate::filter::inputs_should_be_output;
use crate::index::{IndexHeader, load_minimizer_hashes};
use crate::server_common::{FilterRequest, FilterResponse};
use axum::{
    Json, Router,
    extract::DefaultBodyLimit,
    routing::{get, post},
};
use rustc_hash::FxHashSet;

/// Shared index file between endpoint calls.
/// Annoyingly, we have to use an Option as the default/empty FxHashSet is not static
static INDEX: Mutex<Option<FxHashSet<u64>>> = Mutex::new(None);

/// Shared index header between endpoint calls.
/// Initalised to a dummy value, which will be replaced when the index is loaded.
static INDEX_HEADER: Mutex<IndexHeader> = Mutex::new(IndexHeader {
    format_version: 0,
    kmer_length: 0,
    window_size: 0,
});

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

    // build our application with a route
    let app = Router::new()
        // `GET /` goes to `root`
        .route("/", get(root))
        // `GET /index_header` returns the index header
        .route("/index_header", get(index_header))
        // `POST /filter` goes to `filter`
        .route("/should_output", post(should_output))
        // Increase the body limit to 2GB to ensure we don't error on large payloads
        .layer(DefaultBodyLimit::max(2147483648));

    // run our app with hyper, listening globally
    let listener = tokio::net::TcpListener::bind("0.0.0.0:".to_owned() + &port.to_string())
        .await
        .unwrap();
    axum::serve(listener, app).await.unwrap();
}

/// Load the index from the specified path.
fn load_index(index_path: PathBuf) {
    let result = load_minimizer_hashes(&Some(&index_path), &None);
    match result {
        Ok((minimizers, header)) => {
            *INDEX.lock().unwrap() = minimizers;
            *INDEX_HEADER.lock().unwrap() = header;
        }
        Err(e) => {
            eprintln!("Failed to load index: {e}");
            std::process::exit(1);
        }
    }
}

/// Basic root, returing a message indicating the index is loaded
/// Endpoint is `/`
pub async fn root() -> String {
    let index = INDEX.lock();
    match index {
        Ok(index) => {
            let index = index.as_ref().expect("Index not loaded");
            let header = INDEX_HEADER.lock().unwrap();
            format!(
                "Index loaded with {} minimizers and header: {:?}",
                index.len(),
                header
            )
        }
        Err(e) => format!("Error accessing index: {e}"),
    }
}

/// Endpoint to return the header of the loaded index
/// Endpoint is `/index_header`
pub async fn index_header() -> Json<IndexHeader> {
    let header = INDEX_HEADER.lock().unwrap();
    Json(header.clone())
}

/// Endpoint which takes a set of hashes, returning whether they match the index
/// Endpoint is `/should_output`
pub async fn should_output(Json(request): Json<FilterRequest>) -> Json<FilterResponse> {
    let index = INDEX.lock();
    match index {
        Ok(index) => {
            let index = index.as_ref().expect("Index not loaded");
            Json(FilterResponse {
                should_output: inputs_should_be_output(
                    index,
                    &request.input,
                    &request.match_threshold,
                    request.deplete,
                ),
            })
        }
        Err(e) => panic!("Error accessing index: {e}"),
    }
}
