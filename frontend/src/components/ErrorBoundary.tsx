import React, { ErrorInfo } from "react";
import { AppState, SetAppState } from "../appState";

interface Props {
  appState: AppState;
  setAppState: SetAppState;
  children?: React.ReactNode;
}

export class ErrorBoundary extends React.Component<Props> {
  constructor(props: Props) {
    super(props);
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

    if (error) {
      return (
        <div className="app-main centered error-screen">
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
            <pre>{this.props.appState.client.errorInfo?.componentStack}</pre>
          </details>
          <details>
            <summary>App State</summary>
            <pre>{JSON.stringify(this.props.appState, null, 2)}</pre>
          </details>
        </div>
      );
    }

    return this.props.children;
  }
}
