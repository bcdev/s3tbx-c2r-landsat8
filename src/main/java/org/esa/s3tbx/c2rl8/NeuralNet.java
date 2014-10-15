package org.esa.s3tbx.c2rl8;

import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.nn.NNffbpAlphaTabFast;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * @author Marco Peters
 */
public class NeuralNet {


    private ThreadLocal<NNffbpAlphaTabFast> threadLocalNet;

    public static NeuralNet read(InputStream stream) {
       return new NeuralNet(readNeuralNetString(stream));
    }

    private NeuralNet(String nnString) {
        threadLocalNet = new NNffbpAlphaTabFastThreadLocal(nnString);
    }

    public NNffbpAlphaTabFast get() {
        return threadLocalNet.get();
    }

    private static String readNeuralNetString(InputStream neuralNetStream) {
        BufferedReader reader = new BufferedReader(new InputStreamReader(neuralNetStream));
        try {
            String line = reader.readLine();
            final StringBuilder sb = new StringBuilder();
            while (line != null) {
                // have to append line terminator, cause it's not included in line
                sb.append(line).append('\n');
                line = reader.readLine();
            }
            return sb.toString();
        } catch (IOException ioe) {
            throw new OperatorException("Could not initialize neural net", ioe);
        } finally {
            try {
                reader.close();
            } catch (IOException ignore) {
            }
        }
    }

    private static class NNffbpAlphaTabFastThreadLocal extends ThreadLocal<NNffbpAlphaTabFast> {

        private final String nnString;

        public NNffbpAlphaTabFastThreadLocal(String nnString) {
            this.nnString = nnString;
        }

        @Override
        protected NNffbpAlphaTabFast initialValue() {
            try {
                return new NNffbpAlphaTabFast(this.nnString);
            } catch (IOException e) {
                throw new OperatorException("Not able to init neural net", e);
            }
        }
    }
}
