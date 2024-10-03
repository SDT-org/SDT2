import React, { ErrorInfo } from "react";
import { AppState, SetAppState, useSetAppError } from "../appState";
import { Dialog, Modal } from "react-aria-components";

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

    (window as any).lastError = {
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
          error: Error(e.reason),
        },
      };
    });
  }

  getMailTo(error?: Error) {
    let mailTo = `mailto:TODO@TODO?subject=SDT2%20crash`;

    if (error) {
      mailTo = `${mailTo}:%20${error.message}&body=${encodeURI(
        error.stack || "",
      )}`;
    }

    return mailTo;
  }

  override render() {
    const error = this.props.appState.client.error;
    const resetAppError = () =>
      this.props.setAppState((previous) => {
        return {
          ...previous,
          client: {
            ...previous.client,
            error: null,
          },
        };
      });

    if (error) {
      return (
        <Modal isDismissable isOpen={true} onOpenChange={resetAppError}>
          <Dialog>
            <div className="error-modal">
              <h1>Something went wrong.</h1>
              <p>
                Please{" "}
                <a href={this.getMailTo(error)} target="_blank">
                  send us an email
                </a>{" "}
                with the error details.
              </p>
              <details open={true}>
                <summary>{error?.message.toString()}</summary>
                <pre>
                  {this.props.appState.client.errorInfo?.componentStack}
                </pre>
              </details>
              <details>
                <summary>App State</summary>
                <pre>{JSON.stringify(this.props.appState, null, 2)}</pre>
              </details>
              <div className="footer">
                <button onClick={resetAppError}>Close</button>
              </div>
            </div>
          </Dialog>
        </Modal>
      );
    }

    return this.props.children;
  }
}
