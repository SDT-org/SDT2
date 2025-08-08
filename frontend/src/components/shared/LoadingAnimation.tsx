import React from "react";

export const LoadingAnimation = (
  props: React.SVGProps<SVGSVGElement> & { squareSize?: number },
) =>
  React.useMemo(() => {
    const duration = props.dur || "6s";
    const squareSize = props.squareSize || "9";
    const values = props.values || "black; white; black";

    return (
      <svg
        xmlns="http://www.w3.org/2000/svg"
        viewBox="0 0 30 30"
        aria-hidden={true}
        {...props}
      >
        <rect id="one" x="0" y="0" width="9" height="9" fill="black">
          <animate
            attributeName="fill"
            values={values}
            dur={duration}
            repeatCount="indefinite"
            begin="0s"
          />
        </rect>
        <rect x="0" y="10" width="9" height="9" fill="black">
          <animate
            attributeName="fill"
            values={values}
            dur={duration}
            repeatCount="indefinite"
            begin="0.5s"
          />
        </rect>
        <rect x="10" y="10" width="9" height="9" fill="black">
          <animate
            attributeName="fill"
            values={values}
            dur={duration}
            repeatCount="indefinite"
            begin="1s"
          />
        </rect>
        <rect x="0" y="20" width={squareSize} height={squareSize} fill="black">
          <animate
            attributeName="fill"
            values={values}
            dur={duration}
            repeatCount="indefinite"
            begin="1.5s"
          />
        </rect>
        <rect x="10" y="20" width={squareSize} height={squareSize} fill="black">
          <animate
            attributeName="fill"
            values={values}
            dur={duration}
            repeatCount="indefinite"
            begin="2s"
          />
        </rect>
        <rect x="20" y="20" width={squareSize} height={squareSize} fill="black">
          <animate
            attributeName="fill"
            values={values}
            dur={duration}
            repeatCount="indefinite"
            begin="2.5s"
          />
        </rect>
      </svg>
    );
  }, [props]);
