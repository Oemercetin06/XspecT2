from WebApp import app

if __name__ == '__main__':
    
    # Lese die gespeicherten Optionen aus der Datei (oder Datenbank) ein
    try:
        with open('saved_options.txt', 'r') as file:
            saved_options = [line.strip() for line in file]
    except FileNotFoundError:
        pass
    
    app.run(host= '0.0.0.0', port=8000, debug=True, threaded=True)

      # Threading makes Multiple Users at once possible
      # host= '0.0.0.0'
