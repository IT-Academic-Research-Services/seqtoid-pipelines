// src/utils/stream.rs
use std::fs::File;
use std::io::{self, Write};
use tokio::sync::mpsc;
use tokio::sync::broadcast;
use tokio_stream::StreamExt;
use tokio_stream::wrappers::BroadcastStream;


/// Generates any number of output streams from a single input stream.
/// The streams are all asynchronous.
///
/// # Arguments
///
/// * `input_rx' - Receiver stream: tokio::mpsc
/// * 'num_streams' - Number of output streams to generate.
/// 
/// # Returns
/// Vector of output streams.
///
pub async fn t_junction<T>(
    input_rx: mpsc::Receiver<T>,
    num_streams: usize,
) -> Vec<BroadcastStream<T>>
where
    T: Clone + Send + 'static,
{
    let (tx, _rx) = broadcast::channel(100);

    // Create N broadcast receivers first
    let mut streams = Vec::new();
    for _ in 0..num_streams {
        let rx = tx.subscribe();
        streams.push(BroadcastStream::new(rx));
    }

    // Spawn the task after creating receivers
    tokio::spawn(async move {
        let mut rx = input_rx;
        while let Some(item) = rx.recv().await {
            if tx.send(item).is_err() {
                break; // No active receivers
            }
        }
    });

    streams
}


/// Takes a stream 
/// The streams are all asynchronous.
///
/// # Arguments
///
/// * `input_rx' - Receiver stream: tokio::mpsc
/// * 'num_streams' - Number of output streams to generate.
///
/// # Returns
/// Vector of output streams.
///
pub async fn tee_file<T> (input_rx: mpsc::Receiver<T>, out_filename: String) -> BroadcastStream<T>
where
    T: Clone + Send + 'static,
{
    
    let out_streams = t_junction(input_rx, 2).await;

    let mut streams = out_streams.into_iter(); 
    let return_stream = streams.next().unwrap(); 
    let file_stream = streams.next().unwrap();

    
    return_stream 
}
