#![allow(dead_code)]
#![allow(unused_imports)]

use askama::Template;
use futures_util::Future;
use once_cell::sync::Lazy;
use reqwest::Client;
use tempfile::tempdir;

use axum::{
    body::Body,
    extract::{Json, Path, State},
    http::{header, HeaderMap, HeaderName, HeaderValue, StatusCode},
    response::{Html, IntoResponse},
    routing::get,
    routing::post,
    Router,
};
use serde::{de::IntoDeserializer, Deserialize, Serialize};
//use serde::{Deserialize, Serialize};
use tokio::{
    fs::File,
    io::{AsyncBufReadExt, AsyncWriteExt},
    sync::{mpsc::Sender, RwLock},
};

use serde_json::Value;

use tower_sessions::Session;
use tracing::debug;
use tron_app::{
    tron_components::{
        self, tn_future, TnActionExecutionMethod, TnFutureString, TnHtmlResponse, TnSelect,
    },
    AppData,
};
use tron_components::{
    text::TnTextInput, TnButton, TnComponentBaseTrait, TnComponentState, TnComponentValue,
    TnContext, TnContextBase, TnEvent, TnTextArea,
};
//use std::sync::Mutex;
use std::{collections::HashMap, env, fs, io::Read, pin::Pin, str::FromStr, sync::Arc};

use std::collections::HashSet;
use std::process::Command;

static SUPPORTED_GENES: Lazy<HashSet<&'static str>> = Lazy::new(|| {
    let mut set = HashSet::new();
    let supported_genes = [
        "CDC42EP5", "LILRA1", "LILRB1", "LILRB4", "LAIR1", "LENG8", "LILRA2", "LILRA5", "LILRB2",
        "LILRB5", "LAIR2", "LENG9", "LILRA3", "LILRA6", "LILRB3", "TTYH1", "HLA-A", "HLA-B", "HLA-C",
    ];
    supported_genes.into_iter().for_each(|k| {
        set.insert(k);
    });

    set
});

// This is the main entry point of the application
// It sets up the application configuration and state
// and then starts the application by calling tron_app::run
#[tokio::main]
async fn main() {
    let api_routes = Router::<Arc<AppData>>::new().route("/dna2gfe", post(dna2gfe));

    let app_config = tron_app::AppConfigure {
        address: [0, 0, 0, 0],
        log_level: Some("server=debug,tower_http=debug,tron_app=info"),
        http_only: true,
        api_router: Some(api_routes),
        cognito_login: false,
        ..Default::default()
    };
    // set app state
    let app_share_data = tron_app::AppData {
        head: None,
        html_attributes: None,
        context_store: RwLock::new(HashMap::default()),
        session_expiry: RwLock::new(HashMap::default()),
        build_context: Arc::new(Box::new(build_context)),
        build_layout: Arc::new(Box::new(layout)),
    };
    tron_app::run(app_share_data, app_config).await
}

// These functions are used to build the application context,
// layout, and event actions respectively
fn build_context() -> TnContext {
    let context = TnContextBase::default();

    let context = Arc::new(RwLock::new(context));
    TnContext { base: context }
}

#[derive(Template)] // this will generate the code...
#[template(path = "./app_page.html", escape = "none")] // using the template in this path, relative                                    // to the `templates` dir in the crate root
struct AppPageTemplate {}

fn layout(_context: TnContext) -> TnFutureString {
    tn_future! {

        let html = AppPageTemplate {

        };
        html.render().unwrap()
    }
}

#[derive(Serialize, Deserialize, Debug)]
struct DnaSeqInput {
    gene_name: String,
    seq_name: String,
    sequence: String,
}

#[derive(Serialize, Deserialize, Debug)]
struct B12FeatureReq {
    locus: String,
    term: String,
    rank: u32,
    sequence: String,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct B12Response {
    locus: String,
    accession: u32,
    term: String,
    rank: u32,
    sequence: String,
}

#[derive(Serialize, Deserialize, Debug)]
struct Dna2GfeResponse {
    valid_input: bool,
    seq_name: String,
    data: Vec<B12Response>,
    gfe: String,
}

fn split_by_case_transition(input: &str) -> Vec<String> {
    let mut result = Vec::new();
    let mut current_part = String::new();

    for (i, ch) in input.chars().enumerate() {
        // If the character is the first one or no transition has been found
        if i == 0
            || (ch.is_uppercase() && current_part.chars().last().unwrap().is_lowercase())
            || (ch.is_lowercase() && current_part.chars().last().unwrap().is_uppercase())
        {
            // Push the previous part (if exists) and reset current_part
            if !current_part.is_empty() {
                result.push(current_part);
            }
            current_part = String::new();
        }
        current_part.push(ch);
    }

    // Push the last part if exists
    if !current_part.is_empty() {
        result.push(current_part);
    }

    result
}

async fn dna2gfe(
    State(_app_data): State<Arc<AppData>>,
    _session: Session,
    Json(payload): Json<Value>,
) -> Json<Dna2GfeResponse> {
    // let _session_id = if let Some(session_id) = session.id() {
    //     session_id
    // } else {
    //     return (StatusCode::FORBIDDEN, response_headers, Body::default());
    // };

    let dna_seq: DnaSeqInput = serde_json::from_value(payload).unwrap();

    if !SUPPORTED_GENES.contains(&dna_seq.gene_name[..]) {
        return Json(Dna2GfeResponse {
            valid_input: false,
            seq_name: dna_seq.seq_name.clone(),
            data: vec![],
            gfe: "".to_string()
        });
    }

    tracing::info!(target: "tron_app", "dna_seq input:{:?}", dna_seq);

    let work_tmp_dir = tempdir().unwrap();
    tracing::info!(target: "tron_app", "Temporary directory: {:?}", work_tmp_dir);

    let fasta_path = work_tmp_dir.path().join("input.fa");
    tracing::info!(target: "tron_app", "fasta_path: {:?}", fasta_path);
    let mut fasta_file = File::create(fasta_path.clone())
        .await
        .expect("create fasta file error");

    fasta_file
        .write_all(format!(">{}\n", dna_seq.gene_name).as_bytes())
        .await
        .expect("writing fasta file error 0");

    fasta_file
        .write_all(format!("{}\n", dna_seq.sequence).as_bytes())
        .await
        .expect("writing fasta file error 1");

    fasta_file
        .flush()
        .await
        .expect("writing fasta file error 2");

    let miniprot_output = Command::new("/opt/bin/miniprot")
        .args([
            &format!("{}", fasta_path.display()),
            &format!("/opt/ref_data/{}_prot.fa", dna_seq.gene_name),
            "-j",
            "2",
            "--trans",
            "--aln",
            "--max-intron-out",
            "20000",
            "-G",
            "20000",
            "--outs=0.975",
            "--outc=0.8",
            "--gff",
        ])
        .output()
        .expect("failed to execute process");

    tracing::info!(target: "tron_app", "miniprot output:{:?}", miniprot_output);
    //tracing::info!(target: "tron_app", "gff_file: {}", String::from_utf8_lossy(&miniprot_output.stdout[..]));

    let mut lines = miniprot_output.stdout.lines();
    let mut data: Vec<(u32, String, String)> = vec![];
    let mut exon_rank = 1_u32;
    let mut intron_rank = 1_u32;
    while let Some(line) = lines.next_line().await.expect("can't read GFF output") {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 2 {
            continue;
        };
        if parts[0] == "##ATN" {
            let seq = parts[1];
            let seq = seq.replace("-", "");
            let exon_introns = split_by_case_transition(&seq);
            tracing::info!(target: "tron_app", "split out: {:?}", exon_introns);
            exon_introns.into_iter().for_each(|s| {
                if s.starts_with(&['A', 'G', 'C', 'T'][..]) {
                    data.push((exon_rank, "exon".to_string(), s));
                    exon_rank += 1;
                } else {
                    data.push((intron_rank, "intron".to_string(), s));
                    intron_rank += 1;
                }
            });
            break;
        }
    }
    tracing::info!(target: "tron_app", "data = {:?}", data);

    if data.is_empty() {
        return Json(Dna2GfeResponse {
            valid_input: false,
            seq_name: dna_seq.seq_name.clone(),
            data: vec![],
            gfe: "".to_string()
        });
    }

    let req_client = Client::new();

    let mut dna2gfe_out = Dna2GfeResponse {
        valid_input: true,
        seq_name: dna_seq.seq_name.clone(),
        data: vec![],
        gfe: "".into()
    };

    for (rank, term, s) in data.iter() {
        let feature = B12FeatureReq {
            locus: dna_seq.gene_name.clone(),
            term: term.clone(),
            rank: *rank,
            sequence: s.to_uppercase(),
        };
        let response = req_client
            .post("https://feature.b12x.org:443/features") // Replace with your target URL
            .json(&feature) // Serialize the payload to JSON
            .send()
            .await
            .expect("http request error");
        if response.status().is_success() {
            // Deserialize the JSON response
            let response_body: B12Response = response.json().await.expect("b12x request fail");
            // println!("Response: {:?}", response_body);
            dna2gfe_out.data.push(response_body.clone());
            tracing::info!(target: "tron_app", "b12x exon request response {:?}", response_body);
        } else {
            tracing::info!(target: "tron_app", "Failed to send request: {}", response.status());
        }
    }

    let gfe = dna2gfe_out
        .data
        .iter()
        .map(|v| format!("{}", v.accession))
        .collect::<Vec<String>>()
        .join("-");

    let gfe = format!("{}w{}", dna_seq.gene_name, gfe);
    dna2gfe_out.gfe = gfe;

    Json(dna2gfe_out)
}
