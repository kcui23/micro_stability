import { app, shell, BrowserWindow, ipcMain, screen } from 'electron'
import { join } from 'path'
import { electronApp, optimizer, is } from '@electron-toolkit/utils'
import icon from '../../resources/icon.png?asset'
import path from 'path'
import execa from 'execa';
import { chmod } from 'fs/promises';

const rPath = 'r-mac';
// const rpath = path.join(app.getAppPath(), rPath);
const rpath = path.join(process.resourcesPath, rPath);
const libPath = path.join(rpath, 'library');
const rscript = path.join(rpath, 'bin', 'R');
// const rScriptPath = path.join(app.getAppPath(), 'r-scripts');
const rScriptPath = path.join(process.resourcesPath, 'r-scripts');
let rProcess = null;

async function makeExecutable(filePath) {
  if (process.platform === 'darwin' || process.platform === 'linux') {
    try {
      await chmod(filePath, '755');
    } catch (error) {
      console.error('Error making file executable:', error);
    }
  }
}

app.commandLine.appendSwitch('lang', 'en-US');

function createWindow() {
  const { width, height } = screen.getPrimaryDisplay().workAreaSize
  // Create the browser window.
  const mainWindow = new BrowserWindow({
    width: width,
    height: height,
    show: false,
    title: 'Micro-stability',
    autoHideMenuBar: true,
    ...(process.platform === 'linux' ? { icon } : {}),
    webPreferences: {
      preload: join(__dirname, '../preload/index.js'),
      sandbox: false
    }
  })

  mainWindow.on('ready-to-show', () => {
    mainWindow.setTitle('Micro-stability');
    mainWindow.show();
  })
  
  // This handler is fired when the window content attempts to open a new window. 
  // It will open the new window's URL using the default browser and prevent the default behavior 
  // (i.e. opening a new window within the Electron app).
  mainWindow.webContents.setWindowOpenHandler((details) => {
    shell.openExternal(details.url)
    return { action: 'deny' }
  })

  // HMR for renderer base on electron-vite cli.
  // Load the remote URL for development or the local html file for production.
  if (is.dev && process.env['ELECTRON_RENDERER_URL']) {
    mainWindow.loadURL(process.env['ELECTRON_RENDERER_URL'])
  } else {
    mainWindow.loadFile(join(__dirname, '../renderer/index.html'))
  }

  mainWindow.on('close', () => {
    mainWindow.webContents.executeJavaScript(`
      localStorage.removeItem('hasRefreshed');
    `);
  });
}

app.on('ready', () => {
  rProcess = execa(rscript,
    // ['--vanilla', '-f', path.join(app.getAppPath(), 'run-api.R')], { 
    ['--vanilla', '-f', path.join(process.resourcesPath, 'run-api.R')], { 
      env: { // all these vars are sent to the R process as environment variables
             // can be accessed in R with Sys.getenv("VAR_NAME")
        'WITHIN_ELECTRON': '1', 
        'RHOME': rpath,
        'R_HOME_DIR': rpath,
        'R_LIBS': libPath,
        // 'RE_SHINY_PORT': shinyPort,
        // 'RE_SHINY_PATH': shinyAppPath,
        'R_LIBS_USER': libPath,
        'R_LIBS_SITE': libPath,
        'R_LIB_PATHS': libPath,
        'RESOURCES_PATH': process.resourcesPath
      }}).catch((e) => {
          console.error(e)
        })
})

// This method will be called when Electron has finished
// initialization and is ready to create browser windows.
// Some APIs can only be used after this event occurs.
app.whenReady().then(async () => {
  // Set app user model id for windows
  electronApp.setAppUserModelId('com.electron')

  BrowserWindow.getAllWindows().forEach(window => {
    window.webContents.executeJavaScript(`
      localStorage.setItem('hasRefreshed', 'false');
    `);
  });

  // Default open or close DevTools by F12 in development
  app.on('browser-window-created', (_, window) => {
    optimizer.watchWindowShortcuts(window)
  })

  // IPC test
  ipcMain.on('ping', () => console.log('pong'))

  await makeExecutable(rscript);
  createWindow()

  app.on('activate', function () {
    // On macOS it's common to re-create a window in the app when the
    // dock icon is clicked and there are no other windows open.
    if (BrowserWindow.getAllWindows().length === 0) createWindow()
  })
})

// Quit when all windows are closed, except on macOS. There, it's common
// for applications and their menu bar to stay active until the user quits
// explicitly with Cmd + Q.
app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') {
    app.quit()
  }
  try {
    rProcess.kill();
  } catch (e) {
    console.error(e);
  }
})

// In this file you can include the rest of your app"s specific main process
// code. You can also put them in separate files and require them here.
