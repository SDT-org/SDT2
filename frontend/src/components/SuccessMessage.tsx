import { Dialog, Heading, Modal, ModalOverlay } from "react-aria-components";

const SuccessMessage = () => (
  <div className="app-message">
    <svg
      className="app-success"
      xmlns="http://www.w3.org/2000/svg"
      viewBox="0 0 24 24"
      aria-hidden="true"
    >
      <g
        style={{
          fill: "none",
          stroke: "currentcolor",
          strokeWidth: 2,
          strokeLinecap: "round",
          strokeLinejoin: "round",
          strokeMiterlimit: 10,
        }}
      >
        <path d="m8.5 13 3 2 4-6" />
        <path d="M18 22H6a4 4 0 0 1-4-4V6a4 4 0 0 1 4-4h12a4 4 0 0 1 4 4v12a4 4 0 0 1-4 4z" />
      </g>
    </svg>
    <Heading slot="title">Export complete</Heading>
  </div>
);

export const SuccessModal = ({ onClose }: { onClose: () => void }) => (
  <ModalOverlay
    className={"react-aria-ModalOverlay modal-overlay"}
    isOpen
    onOpenChange={onClose}
  >
    <Modal
      isOpen
      onOpenChange={onClose}
      isDismissable={true}
      className={"react-aria-Modal"}
    >
      <Dialog>
        <SuccessMessage />
      </Dialog>
    </Modal>
  </ModalOverlay>
);
