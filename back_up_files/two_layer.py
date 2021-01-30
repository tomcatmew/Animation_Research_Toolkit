import numpy as np
#import matplotlib.pyplot as plt



class TwoLayerNet(object):

    def __init__(self, input_size, hidden_size, output_size, std=1e-4):

        self.params = {}
        self.params['W1'] = std * np.random.randn(input_size, hidden_size)
        self.params['b1'] = np.zeros(hidden_size)
        self.params['W2'] = std * np.random.randn(hidden_size, output_size)
        self.params['b2'] = np.zeros(output_size)
    def loss(self, X, y=None, reg=0.0):

        # Unpack variables from the params dictionary
        W1, b1 = self.params['W1'], self.params['b1']
        W2, b2 = self.params['W2'], self.params['b2']
        N, D = X.shape
        # Compute the forward pass
        scores = None

        # First layer pre-activation
        z1 = X.dot(W1) + b1
        # First layer activation
        a1 = np.maximum(0, z1)
        # Second layer pre-activation
        z2 = a1.dot(W2) + b2
        scores = z2

        # If the targets are not given then jump out, we're done
        if y is None:
            return scores
        # Compute the loss
        loss = None


        mse = ((y - scores) ** 2).mean(axis=0)


        # Backward pass: compute gradients
        grads = {}

        #  my back propagation here
    def train(self, X, y, X_val, y_val,
              learning_rate=1e-3, learning_rate_decay=0.95,
              reg=1e-5, num_iters=100,
              batch_size=200, verbose=False):

        num_train = X.shape[0]
        iterations_per_epoch = max(num_train / batch_size, 1)
        # Use SGD to optimize the parameters in self.model
        loss_history = []
        train_acc_history = []
        val_acc_history = []
        for it in range(num_iters):
            X_batch = None
            y_batch = None

            sample_indices = np.random.choice(num_train, batch_size)
            X_batch = X[sample_indices]
            y_batch = y[sample_indices]

            # Compute loss and gradients using the current minibatch
            loss, grads = self.loss(X_batch, y=y_batch, reg=reg)
            loss_history.append(loss)

            self.params['W1'] += -learning_rate * grads['W1']
            self.params['b1'] += -learning_rate * grads['b1']
            self.params['W2'] += -learning_rate * grads['W2']
            self.params['b2'] += -learning_rate * grads['b2']

            if verbose and it % 100 == 0:
                print('iteration %d / %d: loss %f' % (it, num_iters, loss))
            # Every epoch, check train and val accuracy and decay learning rate.
            if it % iterations_per_epoch == 0:
                # Check accuracy
                train_acc = (self.predict(X_batch) == y_batch).mean()
                val_acc = (self.predict(X_val) == y_val).mean()
                train_acc_history.append(train_acc)
                val_acc_history.append(val_acc)
                # Decay learning rate
                learning_rate *= learning_rate_decay
        return {
            'loss_history': loss_history,
            'train_acc_history': train_acc_history,
            'val_acc_history': val_acc_history,
        }
    def predict(self, X):

        y_pred = None

        z1 = X.dot(self.params['W1']) + self.params['b1']
        a1 = np.maximum(0, z1)  # pass through ReLU activation function
        scores = a1.dot(self.params['W2']) + self.params['b2']
        y_pred = np.argmax(scores, axis=1)

        return y_pred