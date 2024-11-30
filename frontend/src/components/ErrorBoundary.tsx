import React, { type ErrorInfo } from "react";
import { Button, Dialog, Modal } from "react-aria-components";
import { type AppState, type SetAppState, initialAppState } from "../appState";
import { formatBytes } from "../helpers";

interface Props {
  appState: AppState;
  setAppState: SetAppState;
  children?: React.ReactNode;
}

export class ErrorBoundary extends React.Component<Props> {
  constructor(props: Props) {
    super(props);
    this.handlePromiseRejection = this.handlePromiseRejection.bind(this);
  }

  override componentDidMount(): void {
    window.addEventListener("unhandledrejection", this.handlePromiseRejection);
  }

  override componentWillUnmount(): void {
    window.removeEventListener(
      "unhandledrejection",
      this.handlePromiseRejection,
    );
  }

  override componentDidCatch(error: Error, errorInfo: ErrorInfo) {
    this.props.setAppState((previous) => {
      return {
        ...previous,
        client: {
          ...previous.client,
          error,
          errorInfo,
        },
      };
    });

    window.LAST_ERROR = {
      error,
      errorInfo,
    };
  }

  handlePromiseRejection(e: PromiseRejectionEvent) {
    this.props.setAppState((previous) => {
      return {
        ...previous,
        client: {
          ...previous.client,
          error: new Error(e.reason),
          errorInfo: e.reason,
        },
      };
    });
  }

  getMailTo(error?: Error, details = "") {
    let mailTo = "mailto:SDT_admin@proton.me?subject=SDT2%20Issue";

    if (error) {
      mailTo = `${mailTo}:%20${error.message}&body=${encodeURI(details)}`;
    }

    return mailTo;
  }

  getIssueUrl(error?: Error, details = "") {
    const url = new URL("https://github.com/SDT-org/SDT2/issues/new");

    if (!error) {
      return url.toString();
    }

    url.searchParams.set("labels", "bug");
    url.searchParams.set("title", error.message);
    url.searchParams.set("body", details);

    return url.toString();
  }

  override render() {
    const error = this.props.appState.client.error;
    const errorInfo = this.props.appState.client.errorInfo;
    const platform = this.props.appState.platform;
    const stats = this.props.appState.compute_stats;
    const objectToHuman = (obj?: unknown) =>
      Object.entries(obj ?? {})
        .map((v) => `${v[0].toUpperCase()}: ${v[1]}`)
        .join("\n");

    const errorDetails = [
      errorInfo?.stack,
      error?.stack,
      errorInfo?.componentStack,
      objectToHuman({ ...platform, memory: formatBytes(platform.memory) }),
      objectToHuman({
        ...stats,
        available_memory: formatBytes(stats?.available_memory || 0),
        required_memory: formatBytes(stats?.required_memory || 0),
      }),
    ]
      .filter(Boolean)
      .join("\n\n---\n\n");

    const resetAppError = () =>
      this.props.setAppState((previous) => {
        return {
          ...previous,
          client: {
            ...previous.client,
            error: null,
            errorInfo: null,
          },
        };
      });

    const resetApp = () => this.props.setAppState(initialAppState);

    if (error) {
      return (
        <Modal isDismissable isOpen={true} onOpenChange={resetAppError}>
          <Dialog>
            <div className="error-modal">
              <h1>Something went wrong.</h1>
              <p>
                Please{" "}
                <a
                  href={this.getMailTo(error, errorDetails)}
                  target="_blank"
                  rel="noreferrer"
                >
                  send us an email
                </a>{" "}
                with the error details, or{" "}
                <a
                  href={this.getIssueUrl(error, errorDetails)}
                  target="_blank"
                  rel="noreferrer"
                >
                  open an issue
                </a>
                .
              </p>
              <details open={true}>
                <summary>{error?.message.toString()}</summary>
                <pre>{errorDetails}</pre>
              </details>
              <details>
                <summary>App State</summary>
                <pre>{JSON.stringify(this.props.appState, null, 2)}</pre>
              </details>
            </div>
            <div className="actions space-between">
              <Button onPress={resetApp}>Reset</Button>
              <Button onPress={resetAppError}>Continue</Button>
            </div>
          </Dialog>
        </Modal>
      );
    }

    return this.props.children;
  }
}
