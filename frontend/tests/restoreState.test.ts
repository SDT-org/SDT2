import { beforeEach, expect, test } from "bun:test";
import {
  type AppState,
  clientStateKey,
  initialAppState,
} from "../src/appState";
import { restoreClientState } from "../src/restoreClientState";

// Adapted from https://stackoverflow.com/a/26177872
const storageMock = () => {
  let storage: { [key: string]: string } = {};

  return {
    setItem: (key: string, value: string) => {
      storage[key] = value || "";
    },
    getItem: (key: string) => (key in storage ? storage[key] || null : null),
    removeItem: (key: string) => {
      delete storage[key];
    },
    get length() {
      return Object.keys(storage).length;
    },
    key: (i: number) => {
      const keys = Object.keys(storage);
      return keys[i] || null;
    },
    clear: () => {
      storage = {};
    },
  };
};

beforeEach(() => {
  if (storageMock().getItem(clientStateKey)) {
    throw new Error("should not be here");
  }
  global.localStorage = storageMock();
});

test("it returns default state when localStorage key is missing", () => {
  const restored = restoreClientState(initialAppState.client);
  expect(restored).toEqual(initialAppState.client);
});

test("it returns default state if invalid state in localStorage", () => {
  localStorage.setItem(clientStateKey, "{ invalid: true");
  const restored = restoreClientState(initialAppState.client);
  expect(restored).toEqual(initialAppState.client);
});

test("doesn't import invalid state", () => {
  const updatedComputeCores = initialAppState.client.compute_cores + 1;
  const state: AppState["client"] & { test: true } = {
    ...initialAppState.client,
    // biome-ignore lint/suspicious/noExplicitAny: test bad values
    dataView: "badvalue" as any,
    compute_cores: updatedComputeCores,
    test: true,
  };
  localStorage.setItem(clientStateKey, JSON.stringify(state));
  const restored = restoreClientState(initialAppState.client);
  expect(restored.compute_cores).toEqual(updatedComputeCores);
  // biome-ignore lint/suspicious/noExplicitAny: test bad values
  expect((restored as any)?.test).toBeUndefined();
});

test("restores valid state", () => {
  const state = {
    ...initialAppState.client,
    client: { ...initialAppState.client, dataView: "" },
  };
  localStorage.setItem(clientStateKey, JSON.stringify(state));
  const restored = restoreClientState(initialAppState.client);
  expect(restored).toEqual(initialAppState.client);
});
