process WATCH_FOR_POD5S {

    cache "$launchDir/pod5_cache"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    val remote_dir

    script:
    """
    pod5_watcher.py
    """

}