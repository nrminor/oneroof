include {CALL_SLACK_ALERT} from "../modules/call_slack_alert"

workflow SLACK_ALERT {
    take:
    alignment

    main:
    CALL_SLACK_ALERT(
    )
} 