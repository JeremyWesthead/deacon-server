//! Functionality to create a server endpoint which can be used to filter based on a pre-loaded index
use std::path::PathBuf;
use std::sync::OnceLock;

use axum::{
    routing::{get, post},
    Json, Router,
};
use reqwest::header;
use rustc_hash::FxHashSet;
use crate::index::{Index, IndexHeader, load_minimizer_hashes};
use crate::filter::{input_matches_index};
use crate::server_common::{FilterRequest, FilterResponse};

static INDEX: OnceLock<FxHashSet<u64>> = OnceLock::new();
static INDEX_HEADER: OnceLock<IndexHeader> = OnceLock::new();

pub async fn run_server(index_path: PathBuf, port: u16) {
    // initialize tracing
    tracing_subscriber::fmt::init();

    // Load the index before starting the server to ensure it's available for requests
    load_index(index_path).await;

    // build our application with a route
    let app = Router::new()
        // `GET /` goes to `root`
        .route("/", get(root))
        // `GET /index_header` returns the index header
        .route("/index_header", get(index_header))
        // `POST /filter` goes to `filter`
        .route("/is_index_match", post(is_index_match));

    // run our app with hyper, listening globally
    let listener = tokio::net::TcpListener::bind("0.0.0.0:".to_owned() + &port.to_string()).await.unwrap();
    axum::serve(listener, app).await.unwrap();
}

/// Load the index from the specified path.
async fn load_index(index_path: PathBuf) {
    let result = load_minimizer_hashes(&Some(index_path), &None).await;
    match result {
        Ok((minimizers, header)) => {
            INDEX.get_or_init(|| minimizers.unwrap());
            INDEX_HEADER.get_or_init(|| header);
        }
        Err(e) => {
            eprintln!("Failed to load index: {}", e);
            std::process::exit(1);
        }
    }
}

// basic handler that responds with a static string
async fn root() -> String {
    let index = INDEX.get().expect("Index not loaded");
    let header = INDEX_HEADER.get().expect("Index header not loaded");

    format!("Index loaded with {} minimizers and header: {:?}", 
            index.len(), header)
}

/// Endpoint to return the header of the loaded index
async fn index_header() -> Json<IndexHeader> {
    let header = INDEX_HEADER.get().expect("Index header not loaded");
    Json(header.clone())
}

// Endpoint which takes a set of hashes, returning whether they match the index
async fn is_index_match(Json(request): Json<FilterRequest>) -> Json<FilterResponse> {
    Json(FilterResponse {
        index_match: input_matches_index(&INDEX.get().expect("Index not loaded"), &request.input, &request.match_threshold),
    })
}


