console.log('renderer.js loaded')

const information = document.getElementById('info')
information.innerText = `This app is using Chrome (v${window.versions.chrome()}), Node.js (v${window.versions.node()}), and Electron (v${window.versions.electron()})`
information.innerText += "\n"
information.innerText += `Response`

information.style.color = 'green'
information.style.fontSize = '24px'
information.style.fontFamily = 'monospace'
information.style.marginTop = '100px'
information.style.fontWeight = 'bold'
information.style.textShadow = '2px 2px 10px rgba(0,0,0,0.5)'
information.style.textDecoration = 'underline'

information.addEventListener('click', () => {
  information.style.color = 'red'
  information.style.textDecoration = 'none'
  information.style.fontSize = '32px'
  information.style.textShadow = '2px 2px 10px rgba(0,0,0,0.5)'
  information.style.fontStyle = 'italic'
  information.style.fontWeight = 'normal'
})

information.addEventListener('dblclick', () => {
  information.style.color = 'blue'
  information.style.textDecoration = 'underline'
  information.style.fontSize = '16px'
  information.style.textShadow = '2px 2px 10px rgba(0,0,0,0.5)'
  information.style.fontStyle = 'normal'
  information.style.fontWeight = 'bold'
})

// add yellow glow effect when mouseover
information.addEventListener('mouseover', () => {
  information.style.textShadow = '2px 2px 10px rgba(255,255,0,0.5)'
})
information.addEventListener('mouseout', () => {
  information.style.textShadow = '2px 2px 10px rgba(0,0,0,0.5)'
})

const func = async ()=> {
  const response = await window.versions.ping()
  information.innerText += `\n${response}`
  console.log(response)
}
func()