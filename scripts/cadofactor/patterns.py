import abc

class Observable(object):
    """ Implements the subject part of the Observer pattern.
    
    An Observer can register with the subject, the subject then can notify 
    the Observer of state changes by calling the Observer's updateObserver()
    method.
    """
    def __init__(self, *args, **kwargs):
        self.__observers = []
        super().__init__(*args, **kwargs)
    
    def subscribeObserver(self, observer):
        if not observer in self.__observers:
            self.__observers.append(observer)
    
    def unsubscribeObserver(self, observer):
        try:
            self.__observers.remove(observer)
        except ValueError:
            pass
    
    def notifyObservers(self, message):
        for observer in self.__observers:
            observer.updateObserver(message)

class Observer(object):
    """ Implements the client part of the Observer pattern. 
    
    This is an abstract base class; the subclass must implement 
    updateObserver. 
    """
    @abc.abstractmethod
    def updateObserver(self, message):
        # Subclasses should implement this
        return
