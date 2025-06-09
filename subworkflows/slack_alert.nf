include {CALL_SLACK_ALERT} from "../modules/call_slack_alert"

workflow SLACK_ALERT {
    // TODO: add coverage TSVs as a channel here with the collect operator or whatever is needed
    take:
    alignment
    consensus
    variants

    main:
    CALL_SLACK_ALERT(
        alignment,
        consensus,
        variants
    )
} 
