# Nextflow DSL2 Training Example

This directory is a small version of the repository's production layout. It is
intended for learning how the entry point, workflow, modules, and configuration
fit together.

## Structure

```text
training-example/
|-- main.nf
|-- nextflow.config
|-- run.sh
|-- bin/shared/
|   `-- extract_csv.nf
|-- examples/
|   `-- messages.csv
|-- modules/
|   `-- bash/
|       |-- create_message.nf
|       |-- uppercase_message.nf
|       |-- summarize_message.nf
|       |-- gather_summaries.nf
|       `-- bin/
|           `-- create_message.sh
|-- nf-test.config
|-- tests/
|   |-- modules/bash/
|   |   |-- create_message.nf.test
|   |   |-- uppercase_message.nf.test
|   |   |-- summarize_message.nf.test
|   |   `-- gather_summaries.nf.test
|   `-- workflows/
|       `-- training.nf.test
`-- workflows/
    `-- training.nf
```

The data flow is:

```text
CSV rows -> (SampleID, Message) tuples -> CREATE_MESSAGE
         -> (SampleID, message file) -> UPPERCASE_MESSAGE
         -> (SampleID, uppercase file) -> SUMMARIZE_MESSAGE
         -> collect all summary files -> GATHER_SUMMARIES
```

The CSV parser validates the `SampleID,Message` headers and converts each row
into a tuple. Every module receives and emits the sample ID alongside its data,
which keeps outputs associated with the correct sample while rows run in
parallel. Each module contains one process and exposes a named output with
`emit`. The workflow imports the parser and modules, then connects their
channels. `main.nf` launches the workflow. `nextflow.config` supplies defaults
and publishes only the final `*.summary.txt` files to `results/`.

`CREATE_MESSAGE` also demonstrates how a process can call a helper script from
its module directory with `${moduleDir}`. The helper script stays next to the
module instead of putting the shell implementation in the process definition.

## Run

Nextflow must be installed or loaded as a module. From this directory, run:

```bash
bash run.sh
```

Or invoke Nextflow directly:

```bash
nextflow run main.nf -profile local --csv_input examples/messages.csv
```

After a successful run, inspect the three files in `results/`. Add another row
to `examples/messages.csv` to see Nextflow schedule another set of tasks. Use
`nextflow clean -f` when the intermediate work files are no longer needed.

## Test

If `nf-test` is installed, run the module and workflow examples with:

```bash
nf-test test
```

The module tests run each process independently. The workflow test runs the
CSV parser and all three modules together, then checks that all three sample
IDs are emitted.

The workflow also demonstrates scatter/gather. Each CSV row scatters through
the first three modules independently. `collect()` waits for all summary files,
then `GATHER_SUMMARIES` runs once and combines them into `results/all_summaries.txt`.