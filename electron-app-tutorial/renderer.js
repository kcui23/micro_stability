console.log('renderer.js loaded')

const information = document.getElementById('info')
information.innerText = `This app is using Chrome (v${window.versions.chrome()}), Node.js (v${window.versions.node()}), and Electron (v${window.versions.electron()})`

information.style.fontSize = '24px'
information.style.fontFamily = 'monospace'
information.style.textShadow = '2px 2px 10px rgba(0,0,0,0.5)'
information.style.textDecoration = 'underline'


// add yellow glow effect when mouseover
information.addEventListener('mouseover', () => {
  information.style.textShadow = '2px 2px 10px rgba(255,255,0,0.5)'
})
information.addEventListener('mouseout', () => {
  information.style.textShadow = '2px 2px 10px rgba(0,0,0,0.5)'
})

const func = async ()=> {
  const response = await window.versions.ping()
  information.innerText += `\nResponse from 'ping': ${response}`
  console.log(response)
}
func()