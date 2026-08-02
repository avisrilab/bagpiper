//! Reusable parallel feeder-consumer.
//!
//! The caller thread produces items via `feed` (so a non-Send source like a FASTQ reader stays on
//! the caller thread); a pool of `workers` threads applies `process`, each holding its own `init`
//! state (e.g. reused buffers); a dedicated writer thread applies `write` to every result. Results
//! arrive in worker-completion order (unordered). Bounded channels give backpressure. Returns the
//! writer's value. Stages that mutate shared output stay single-threaded in the writer, so no locks.

use std::thread;

use crossbeam_channel::{bounded, Receiver};

/// Worker count: the machine's parallelism minus the feeder and writer threads.
pub fn default_workers() -> usize {
    thread::available_parallelism()
        .map(|n| n.get().saturating_sub(2).max(1))
        .unwrap_or(1)
}

/// Run the feeder-consumer. `feed` runs on the calling thread; `process` on `workers` threads, each
/// with its own `init` state; `write` on one writer thread, consuming the result receiver and
/// returning a value that is propagated back to the caller.
pub fn run<I, O, S, R>(
    mut feed: impl FnMut() -> Option<I>,
    workers: usize,
    init: impl Fn() -> S + Sync,
    process: impl Fn(&mut S, I) -> O + Sync,
    write: impl FnOnce(Receiver<O>) -> R + Send,
) -> R
where
    I: Send,
    O: Send,
    R: Send,
{
    let workers = workers.max(1);
    let (item_tx, item_rx) = bounded::<I>(workers * 4);
    let (res_tx, res_rx) = bounded::<O>(workers * 4);

    thread::scope(|s| {
        for _ in 0..workers {
            let item_rx = item_rx.clone();
            let res_tx = res_tx.clone();
            let init = &init;
            let process = &process;
            s.spawn(move || {
                let mut state = init();
                for item in item_rx {
                    if res_tx.send(process(&mut state, item)).is_err() {
                        break;
                    }
                }
            });
        }
        drop(item_rx);
        drop(res_tx);

        let writer = s.spawn(move || write(res_rx));

        while let Some(item) = feed() {
            if item_tx.send(item).is_err() {
                break;
            }
        }
        drop(item_tx);

        writer.join().expect("writer thread panicked")
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn processes_every_item_unordered() {
        let mut n = 0u32;
        let sum: u64 = run(
            || (n < 1000).then(|| {
                n += 1;
                n
            }),
            4,
            || (),
            |_: &mut (), x: u32| (x as u64) * (x as u64),
            |rx: Receiver<u64>| rx.into_iter().sum(),
        );
        let expect: u64 = (1..=1000u64).map(|x| x * x).sum();
        assert_eq!(sum, expect);
    }
}
