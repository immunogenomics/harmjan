package org.anarres.parallelgzip;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * @author shevek
 */
public class ParallelGZIPEnvironment {

	public static ExecutorService getSharedThreadPool() {
		return ThreadPoolHolder.EXECUTOR;
	}

	private static class ThreadPoolHolder {

		private static final ExecutorService EXECUTOR = Executors.newCachedThreadPool();
	}
}
