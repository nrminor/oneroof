process WATCH_FOR_POD5S {

    cache params.pod5_staging

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path config

    output:
    path "*.pod5"

    script:
    """
    pod5_watcher.py --host_config ${config}
    """

}
