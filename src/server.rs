//! Functionality to create a server endpoint which can be used to filter based on a pre-loaded index
use std::path::PathBuf;
use std::sync::OnceLock;

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
static INDEX: OnceLock<FxHashSet<u64>> = OnceLock::new();
/// Shared index header between endpoint calls.
static INDEX_HEADER: OnceLock<IndexHeader> = OnceLock::new();

/// Starts the server with the given index path and port.
/// To log the server's connections, set `RUST_LOG=trace` in your environment variables.
pub async fn run_server(index_path: PathBuf, port: u16) {
    // initialize tracing
    tracing_subscriber::fmt::init();

    // Load the index before starting the server to ensure it's available for requests
    load_index(index_path);

    eprintln!(
        "Index loaded with {} minimizers and header: {:?}",
        INDEX.get().expect("Index not loaded").len(),
        INDEX_HEADER.get().expect("Index header not loaded")
    );

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
            INDEX.get_or_init(|| minimizers.unwrap());
            INDEX_HEADER.get_or_init(|| header);
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
    let index = INDEX.get().expect("Index not loaded");
    let header = INDEX_HEADER.get().expect("Index header not loaded");

    format!(
        "Index loaded with {} minimizers and header: {:?}",
        index.len(),
        header
    )
}

/// Endpoint to return the header of the loaded index
/// Endpoint is `/index_header`
pub async fn index_header() -> Json<IndexHeader> {
    let header = INDEX_HEADER.get().expect("Index header not loaded");
    Json(header.clone())
}

/// Endpoint which takes a set of hashes, returning whether they match the index
/// Endpoint is `/should_output`
pub async fn should_output(Json(request): Json<FilterRequest>) -> Json<FilterResponse> {
    Json(FilterResponse {
        should_output: inputs_should_be_output(
            INDEX.get().expect("Index not loaded"),
            &request.input,
            &request.match_threshold,
            request.deplete,
        ),
    })
}
