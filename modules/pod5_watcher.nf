process WATCH_FOR_POD5S {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    script:
    println("Hi mom!")

}